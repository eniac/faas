lukev@seas.upenn.edu 
Last updated 08/22/2015

Overview:
    This document describes how to set up your Amazon AWS environment for Factoring as a Service.

Steps:
1. Set up and configure AWS CLI using the following instructions:

    http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-set-up.html

    Tips: 
        - Look at spot price history for the instance type that you want to use in different regions, and choose one that has a reasonably stable price history.
        - Choose a region that is physically close to you for reduced latency to the instances.

2. Create VPC (Virtual Private Cloud)
    $ aws ec2 create-vpc --cidr-block 10.0.0.0/20

    2.1. Create Internet Gateway
    2.2 Attach Internet Gateway to VPC
    2.3 Create route table entry for Internet Gateway
    2.4 Make sure DNS hostnames and resolution are enabled

3. Create security group using the id of the previously created VPC.
    $ aws ec2 create-security-group --group-name <my-security-group> --vpc-id <my-vpc-id> --description "my security group"

    Tips:
        - Set the permissions on your security group to allow connections from the ports used by slurm and nfs.
        - Set permissions on your security group to allow ssh on port 22 from the IP address of the machine used to launch the factoring job.

4. Create subnets using id of previously created VPC for each availability zone that you might want to launch nodes into.
    $ aws ec2 create-subnet --vpc-id <my-vpc-id> --cidr-block 10.0.0.0/24 --availability-zone <availability-zone>

    Tips:
        - Execute this command for each availability zone in your region with different subnets, so you can choose your availability zones based on current spot prices.
        - Make sure the subnet assigns a public IP by default
        - Ensure that the route table for the subnet allows access to the outside world.

5. Create placement group
    $ aws ec2 create-placement-group --group-name <my-placement-group> --strategy cluster

6. Create a key pair
    $ aws ec2 create-key-pair --key-name <my-key-pair> >> ~/.ssh/<my-key-pair>.pem
    $ chmod 400 ~/.ssh/<my-key-pair>.pem
