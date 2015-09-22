#!/bin/bash

# Copy params file over to master node
ansible-playbook -i ./ec2.py factor.yml
