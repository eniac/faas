#!/bin/bash

# Terminate if any one command exits with non-zero status
set -e

# Build AMI

## Start ami_builder node
ansible-playbook -i inventory/aws_ec2.yml build-start.yml

## Build image
ansible-playbook -i inventory/aws_ec2.yml build.yml --tags custom --extra-vars="test=yes"

## Create AMI from running instance and then terminate instance
ansible-playbook -i inventory/aws_ec2.yml build-finish.yml 

# Launch cluster

## Launch the master node first so that the NFS and Slurm controller are up before the mpi and slave nodes try to connect.
ansible-playbook -i inventory/aws_ec2.yml launch.yml --tags master --extra-vars="test=yes" 

## Request spot instances in parallel so instances are launched more quickly. 
parallel -u ansible-playbook -i inventory/aws_ec2.yml launch.yml --tags {} --extra-vars="test=yes" ::: mpi slave0 #slave1 slave2 slave3 # <-- uncomment if you're using multiple launch groups

# Start factorization

ansible-playbook -i inventory/aws_ec2.yml factor.yml --extra-vars="test=yes"
