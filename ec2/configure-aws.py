#!/usr/bin/python

import subprocess
import json
import shlex
import os
import sys

def run_command(command):
    args = shlex.split(command)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out,err

print('--- Starting AWS Configuration ---')

# Figure out which region we are in
region = None
aws_config_file = os.path.join(os.path.expanduser('~'), '.aws/config')
with open(aws_config_file, 'r') as f:
    for line in f:
        if line.lower().startswith('region'):
            region = line.split('=')[1].strip()
if region == None:
    print('Unable to find region in ~/.aws/config')
    sys.exit()

print('Using region {} found in ~/.aws/config'.format(region))

ec2_regions = ['us-east-1', 'us-west-1', 'us-west-2', 'eu-west-1', 'eu-central-1', 'ap-southeast-1', 'ap-southeast-2', 'ap-northeast-1', 'sa-east-1']

if region not in ec2_regions:
    print("Unrecognized EC2 region: {}".format(region))
    sys.exit()

# Get subnets in region that support the spot instance types that we will launch,
# since not every availability zone supports all instance types.
from datetime import datetime
now = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
instance_type = 'c4.8xlarge'
history = json.loads(run_command('aws ec2 describe-spot-price-history --instance-types {instance_type} --start-time {now} --end-time {now} --product-description "Linux/UNIX (Amazon VPC)"'.format(now=now, instance_type=instance_type))[0]).get('SpotPriceHistory')
if not history:
    print("Unable to check spot price history for instance type {instance_type} in region {region}".format(instance_type=instance_type, region=region))
    sys.exit()
availability_zones = [x['AvailabilityZone'] for x in history]

# Pick a subnet for each availability zone.
zone_subnet_map = {}
for index,zone in enumerate(sorted(availability_zones)):
    zone_subnet_map[zone] = '10.0.{}.0/24'.format(index)

# Create VPC (http://docs.aws.amazon.com/cli/latest/reference/ec2/create-vpc.html). Note: this may fail if you already have too many VPCs.
cidr_block = '10.0.0.0/16'
vpc_id = json.loads(run_command('aws ec2 create-vpc --cidr-block {cidr_block}'.format(cidr_block=cidr_block))[0])['Vpc']['VpcId']
print("Created VPC {vpc_id}".format(vpc_id=vpc_id))

# Create Internet Gateway (http://docs.aws.amazon.com/cli/latest/reference/ec2/create-internet-gateway.html)
ig_id = json.loads(run_command('aws ec2 create-internet-gateway')[0])['InternetGateway']['InternetGatewayId']
print("Created Internet Gateway {ig_id}".format(ig_id=ig_id))

# Attach Internet Gateway to VPC, after waiting to make sure that the VPC is available
vpc_ready = False
while not vpc_ready:
    if not vpc_ready:
        vpc_ready = (json.loads(run_command('aws ec2 describe-vpcs --vpc-ids {vpc_id}'.format(vpc_id=vpc_id))[0])['Vpcs'][0]['State'] == 'available')
run_command('aws ec2 attach-internet-gateway --internet-gateway-id {ig_id} --vpc-id {vpc_id}'.format(ig_id=ig_id, vpc_id=vpc_id))
print("Attached Internet Gateway {ig_id} to VPC {vpc_id}".format(ig_id=ig_id, vpc_id=vpc_id))

# Enable DNS resolution and hostnames for VPC
run_command('aws ec2 modify-vpc-attribute --vpc-id {vpc_id} --enable-dns-support'.format(vpc_id=vpc_id))
run_command('aws ec2 modify-vpc-attribute --vpc-id {vpc_id} --enable-dns-hostnames'.format(vpc_id=vpc_id))
print("Enabled DNS support and hostnames for VPC")

# Create Security Group
sg_name = "faas"
sg_id = json.loads(run_command('aws ec2 create-security-group --group-name \"{sg_name}\" --vpc-id {vpc_id} --description \"{sg_desc}\"'.format(vpc_id=vpc_id, sg_name=sg_name, sg_desc="faas security group"))[0])['GroupId']
print("Created security group {sg_name} with id {sg_id}".format(sg_name=sg_name, sg_id=sg_id))

# Add ingress rules to security group to allow SSH access. TODO: Modify this to be more restrictive. Ports 22 and 9001 are probably the only ones that need to remain open for SSH and supervisor, and only from the IP address of the control machine.
run_command('aws ec2 authorize-security-group-ingress --group-id {sg_id} --protocol all --port=all --cidr 0.0.0.0/0'.format(sg_id=sg_id))
print("Added ingress rule to security group")

# Create route table for VPC
rt_id = json.loads(run_command('aws ec2 create-route-table --vpc-id {vpc_id}'.format(vpc_id=vpc_id))[0])['RouteTable']['RouteTableId']
# Add entry to route table specifying gateway
run_command('aws ec2 create-route --route-table-id {rt_id} --destination-cidr-block 0.0.0.0/0 --gateway-id {ig_id}'.format(rt_id=rt_id, ig_id=ig_id)) 
print("Created route table {rt_id} with entry to gateway {ig_id}".format(rt_id=rt_id, ig_id=ig_id))

# Create subnets for each availability zone that you want to use. Since we are added all possible ip addresses from these subnets to the Slurm config, use small subnets (/24 or smaller). 
vpc_subnets = []
for az,cidr in zone_subnet_map.items():
    subnet = (json.loads(run_command('aws ec2 create-subnet --vpc-id {vpc_id} --cidr-block {cidr} --availability-zone {az}'.format(vpc_id=vpc_id, cidr=cidr, az=az))[0])['Subnet'])
    # Associate previously created route table with subnet
    run_command('aws ec2 associate-route-table --route-table-id {rt_id} --subnet-id {sn_id}'.format(rt_id=rt_id, sn_id = subnet['SubnetId']))
    vpc_subnets.append(subnet)
    # Auto-assign public ip in subnet
    run_command('aws ec2 modify-subnet-attribute --subnet-id {sn_id} --map-public-ip-on-launch'.format(sn_id=subnet['SubnetId']))
    print('Created subnet {subnet_id} in availability zone {az} with cidr block {cidr}'.format(subnet_id = subnet['SubnetId'], az=subnet['AvailabilityZone'], cidr=subnet['CidrBlock'], rt_id=rt_id))

# Create placement group
placement_group = 'faas-cluster'
try:
    json.loads(run_command('aws ec2 describe-placement-groups --group-name {}'.format(placement_group))[0])
    print('Placement group {} already exists'.format(placement_group))
except Exception as e:
    run_command('aws ec2 create-placement-group --group-name {} --strategy cluster'.format(placement_group))
    print('Created placement group {}'.format(placement_group))

# Create key pair
key_pair = 'faas'
ssh_dir = os.path.dirname('{home}/.ssh/'.format(home=os.path.expanduser('~')))
try: 
    os.makedirs(ssh_dir)
except OSError:
    if not os.path.isdir(ssh_dir):
        raise
os.chmod(ssh_dir, 0o700)
key_file = os.path.join(ssh_dir, '{key_pair}.pem'.format(key_pair=key_pair))
try:
    json.loads(run_command('aws ec2 describe-key-pairs --key-name={}'.format(key_pair))[0])
    print('Key pair {} already exists'.format(key_pair))
    if not os.path.exists(key_file):
        print('Key file {key_file} does not exist locally, but it is present online. Please either remove the key pair {{ key_pair }} via the AWS console, or move the key to the appropriate location locally.'.format(key_file=key_file, key_pair=key_pair))
        sys.exit()

except Exception as e:

    if os.path.exists(key_file):
        print('Key file {key_file} already exists. Please move or rename it.'.format(key_file=key_file))
        sys.exit()

    key_str = json.loads(run_command('aws ec2 create-key-pair --key-name {key_pair}'.format(key_pair=key_pair))[0])['KeyMaterial']
    with open(key_file, 'w') as f:
        f.write(key_str)
    os.chmod(key_file, 0o600)
    print('Created key pair {key_pair} at {key_file}'.format(key_pair=key_pair, key_file=key_file))

# Find an Ubuntu base image to use for this availability region.
base_image = json.loads(run_command('aws ec2 describe-images --filters Name=name,Values=ubuntu/images/hvm-ssd/ubuntu-trusty-14.04-amd64-server-20150325')[0])['Images'][0]['ImageId']

# Write variables to file in yaml format
var_file = 'vars/ec2.yml'
print('Writing EC2 variables to file {var_file}'.format(var_file=var_file))
with open(var_file, 'w') as f:
    f.write('# This file is created and will be overwritten by configure-aws.py.\n')
    f.write('---\n')
    f.write('ec2:\n')
    f.write('    region: {}\n'.format(region))
    f.write('    vpc_id: {}\n'.format(vpc_id))
    f.write('    cidr_block: {}\n'.format(cidr_block))
    f.write('    device_name: /dev/sda1\n')
    f.write('    security_group: {}\n'.format(sg_name))
    f.write('    placement_group: {}\n'.format(placement_group))
    f.write('    ssh_key: {}\n'.format(key_pair))
    f.write('    base_image: {}\n'.format(base_image))
    f.write('    custom_image: \'\'\n')
    f.write('    image_name: faas image\n')
    f.write('    master_private_ip: 10.0.0.4\n')
    f.write('    subnets:\n')
    for sn in sorted(vpc_subnets, key=lambda sn: sn['CidrBlock']):
        f.write('        - {{ availability_zone: {az}, subnet_id: {sn_id}, cidr_block: {cidr} }}\n'.format(az=sn['AvailabilityZone'], sn_id=sn['SubnetId'], cidr=sn['CidrBlock']))

print('--- Finished AWS Configuration ---')
