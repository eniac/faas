#!/bin/bash

# Launch processes in parallel so instances are launched more quickly. We launch the master node first so that NFS and the Slurm controller are up before the MPI and slave nodes try to connect.

# To run the test
#ansible-playbook -i ./ec2.py launch.yml --tags master --extra-vars="test=yes" 
#parallel -u ansible-playbook -i ./ec2.py launch.yml --tags {} --extra-vars="test=yes" ::: master mpi slave0 slave1 slave2 slave3

ansible-playbook -i ./ec2.py launch.yml --tags master
parallel -u ansible-playbook -i ./ec2.py launch.yml --tags {} ::: mpi slave0 slave1 slave2 slave3
