#!/bin/bash

#ansible-playbook -i ./ec2.py factor.yml --extra-vars="test=yes"
ansible-playbook -i ./ec2.py factor.yml
