#!/bin/bash

# Terminate if any one command exits with non-zero status
set -e

ansible-playbook -i inventory/aws_ec2.yml check-cpu-usage.yml
