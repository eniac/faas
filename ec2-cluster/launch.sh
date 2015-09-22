#!/bin/bash

# Launch these processes in the background so that they can all run simultaneously.
ansible-playbook -i ./ec2.py launch-cluster.yml --tags=master &
ansible-playbook -i ./ec2.py launch-cluster.yml --tags=mpi &
ansible-playbook -i ./ec2.py launch-cluster.yml --tags=slave &
wait
