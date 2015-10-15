#!/bin/bash

# Terminate if any one command exits with non-zero status
set -e

# Build AMI

## Start ami_builder node
ansible-playbook -i ./ec2.py build-start.yml

## Build image
ansible-playbook -i ./ec2.py build.yml --tags custom

## Create AMI from running instance and then terminate instance
ansible-playbook -i ./ec2.py build-finish.yml 

# Launch cluster

## Launch the master node first so that the NFS and Slurm controller are up before the mpi and slave nodes try to connect.
ansible-playbook -i ./ec2.py launch.yml --tags master

## Request spot instances in parallel so instances are launched more quickly. 
parallel -u ansible-playbook -i ./ec2.py launch.yml --tags {} ::: mpi slave0 #slave1 slave2 slave3 # <-- uncomment if you're using multiple launch groups

# Start factorization

ansible-playbook -i ./ec2.py factor.yml
