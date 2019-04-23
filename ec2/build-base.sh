#!/bin/bash

# This script takes 20-30 minutes to run, and creates a base image.
set -e
ansible-playbook -i ./ec2.py build-start.yml
ansible-playbook -i ./ec2.py build.yml --tags base
ansible-playbook -i ./ec2.py build-finish.yml  --tags base
