#!/usr/bin/python

# LTV
# Last Edit: 10.07.15
# Create and configure a VPC for FaaS

import subprocess
import json
import shlex
import os

def run_command(command):
    args = shlex.split(command)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out,err

print('--- Starting AWS Configuration ---')

# Create VPC (http://docs.aws.amazon.com/cli/latest/reference/ec2/create-vpc.html)
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
print("Created security group \"{sg_name}\"; {sg_id}".format(sg_name=sg_name, sg_id=sg_id))

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
for az,cidr in [('us-east-1b', '10.0.0.0/24'), ('us-east-1c', '10.0.1.0/24'), ('us-east-1d', '10.0.2.0/24'), ('us-east-1e', '10.0.3.0/24')]:
    subnet = (json.loads(run_command('aws ec2 create-subnet --vpc-id {vpc_id} --cidr-block {cidr} --availability-zone {az}'.format(vpc_id=vpc_id, cidr=cidr, az=az))[0])['Subnet'])
    # Associate previously created route table with subnet
    run_command('aws ec2 associate-route-table --route-table-id {rt_id} --subnet-id {sn_id}'.format(rt_id=rt_id, sn_id = subnet['SubnetId']))
    vpc_subnets.append(subnet)
    # Auto-assign public ip in subnet
    run_command('aws ec2 modify-subnet-attribute --subnet-id {sn_id} --map-public-ip-on-launch'.format(sn_id=subnet['SubnetId']))
    print('Created subnet {subnet_id} in availability zone {az} with cidr block {cidr} with associated route table {rt_id}'.format(subnet_id = subnet['SubnetId'], az=subnet['AvailabilityZone'], cidr=subnet['CidrBlock'], rt_id=rt_id))

# Create placement group
placement_group = 'faas-cluster2'
run_command('aws ec2 create-placement-group --group-name {placement_group} --strategy cluster'.format(placement_group=placement_group))
print('Created placement group {placement_group}'.format(placement_group=placement_group))

# Create key pair
key_pair = 'faas2'
key_str = json.loads(run_command('aws ec2 create-key-pair --key-name {key_pair}'.format(key_pair=key_pair))[0])['KeyMaterial']
key_file = '{home}/.ssh/{key_pair}.pem'.format(home=os.path.expanduser('~'), key_pair=key_pair)
with open(key_file, 'w') as f:
    f.write(key_str)
os.chmod(key_file, 0o600)
print('Created key pair {key_pair} at {key_file}'.format(key_pair=key_pair, key_file=key_file))

# Write variables to file in yaml format
var_file = 'vars/ec2.yml'
print('Writing EC2 variables to file {var_file}'.format(var_file=var_file))
with open(var_file, 'w') as f:
    f.write('# This file is created and will be overwritten by configure-aws.py.\n')
    f.write('---\n')
    f.write('ec2:\n')
    f.write('    region: us-east-1\n')
    f.write('    vpc_id: {vpc_id}\n'.format(vpc_id=vpc_id))
    f.write('    cidr_block: {cidr_block}\n'.format(cidr_block=cidr_block))
    f.write('    device_name: /dev/sda1\n')
    f.write('    security_group: {sg_name}\n'.format(sg_name=sg_name))
    f.write('    placement_group: {placement_group}\n'.format(placement_group=placement_group))
    f.write('    ssh_key: {key_pair}\n'.format(key_pair=key_pair))
    f.write('    base_image: ami-d05e75b8\n') # ami-d05e75b8
    f.write('    custom_image: \'\'\n')
    f.write('    image_name: faas image\n')
    f.write('    master_private_ip: 10.0.0.4\n')
    f.write('    subnets:\n')
    for sn in vpc_subnets:
        f.write('        - {{ availability_zone: {az}, subnet_id: {sn_id}, cidr_block: {cidr} }}\n'.format(az=sn['AvailabilityZone'], sn_id=sn['SubnetId'], cidr=sn['CidrBlock']))

print('--- Finished AWS Configuration ---')
