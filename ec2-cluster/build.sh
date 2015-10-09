#!/bin/bash

# Cause script to terminate if any one command exits with non-zero status
set -e

# Start ami_builder node
ansible-playbook -i ./ec2.py build-start.yml

#ansible-playbook -i ./ec2.py build.yml
#ansible-playbook -i ./ec2.py build.yml --tags custom --extra-vars="test=yes"
ansible-playbook -i ./ec2.py build.yml --tags custom
#ansible-playbook -i ./ec2.py build.yml --tags base

# Create AMI from running instance and then terminate instance
ansible-playbook -i ./ec2.py build-finish.yml 
