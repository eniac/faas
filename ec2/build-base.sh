#!/bin/bash

# This script takes 20-30 minutes to run, and creates a base image.
set -e
ansible-playbook -i inventory/aws_ec2.yml build-start.yml
ansible-playbook -i inventory/aws_ec2.yml build.yml --tags base
ansible-playbook -i inventory/aws_ec2.yml build-finish.yml  --tags base
