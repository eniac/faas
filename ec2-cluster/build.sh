#!/bin/bash

# Cause script to terminate if any one command exits with non-zero status
set -e

ansible-playbook -i ./ec2.py build-start.yml
ansible-playbook -i ./ec2.py build.yml
ansible-playbook -i ./ec2.py build-finish.yml 
