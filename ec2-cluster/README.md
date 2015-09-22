lukev@seas.upenn.edu 
Last updated 08/22/2015

Overview:

    The purpose of these Ansible playbooks are to do the following:

        1. Build a custom AMI for the EC2 environment.
        2. Launch an EC2 cluster configured with slurm and NFS.
        3. Create a custom cado-nfs parameter file and copy to master node.
        4. Start a factoring job monitored using supervisor.

Instructions:

1. Command machine setup:

    "Follow EC2_SETUP.md to set up your EC2 enviroment and AWS CLI."
    
    "Edit ec2 configurations in ec2_variables.py"
        ansible_ssh_private_key_file: <path_to_ec2_key> 
        ec2:
          security_group: <ec2 security group>
          placement_group: <ec2 placement group>
          ssh_key: <ec2 key name>
          aws_region: <region>
          aws_zone: <availability-zone>
          vpc_subnet_id: <subnet-id>
          vpc_id: <vpc-id>
          image_name: <my-image-name>
          ami_description: contains slurm, nfs, cado-nfs, supervisor

    "Install Ansible on command machine (e.g., your workstation)"
        http://docs.ansible.com/ansible/intro_installation.html#installation

2. Build an AMI:

    "build.sh builds a new EC2 image by calling the three playbooks below."
    $ build.sh
    
    "launch a new ami_builder instance"
    $ ansible-playbook -i ./ec2.py build-start.yml 

    "install packages on ami_builder instance"
    $ ansible-playbook -i ./ec2.py build.yml 
    OR, to build an AMI without msieve and cado installed
    $ ansible-playbook -i ./ec2.py build.yml --skip-tags=msieve,cado

    "remove a previous AMI (if found), create a new AMI from the running instance, and terminate the ami_builder instance"
    $ ansible-playbook -i ./ec2.py build-finish.yml

3. Launch a cluster:

    "Edit cluster configurations in ec2_variables.py"
        my_tag: <cluster-tag-prefix>
        master:
            type: <instance-type>
            cores: <instance-cores>
            volume_size: <instance-volume-size>
        mpi:
            type: <instance-type>
            cores: <instance-cores>                         # cores must match the number of virtual CPUs for instance type
            volume_size: <instance-volume-size>
            count: <number-of-instances>
            spot_price: <max-spot-bid>                      # leave commented out to request dedicated instance
            spot_wait_timeout: <spot-bid-timeout>           # leave commented out to request dedicated instance           
        slave:
            type: <instance-type>
            cores: <instance-cores>                         # cores must match the number of virtual CPUs for instance type
            volume_size: <instance-volume-size>
            count: <number-of-instances>
            spot_price: <max-spot-bid>                      # leave commented out to request dedicated instance
            spot_wait_timeout: <spot-bid-timeout>           # leave commented out to request dedicated instance           

    "launch a cluster on ec2"
    $ ansible-playbook -i ./ec2.py launch-cluster.yml

4. Start a factoring job:
    
    "Edit factoring configurations in ec2_variables.py"
        supervisor:
            username: <supervisor-username>
            password: <supervisor-password>
            email: <email-for-crash-reports>                # Ensure that your spam filter doesn't prevent incoming emails from crashmail
        cado:
            name: <task_name>
            N: <number-to-factor>
            <other params>: <values>

    "Copy cado parameters and supervisor configurations to master node"
    $ ansible-playbook -i ./ec2.py factor.yml

    "Start server process via supervisor"
    $ supervisorctl -s http://<master-ip>:<supervisor-port>
        username: <supervisor-username>
        password: <supervisor-password>
        supervisor> start server

    "terminate instances (should only be required if the factoring job doesn't complete as expected)"
    $ ansible-playbook -i ./ec2.py terminate-cluster.yml
