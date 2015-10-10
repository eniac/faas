# FaaS
The purpose of the FaaS (Factoring as a Service) project is to create a tool that allows anyone to factor a 512-bit RSA key with minimal time and cost. To do this, we have developed a framework for launching an compute cluster on Amazon EC2 and running a parallelized implementation of the Number Field Sieve (NFS) algorithm to factor keys. For more information about the project, see [the project webpage](https://www.cis.upenn.edu/~securlab/projects/faas/).

The tool was built using the AWS region 'us-east'. Run in any other region at your own risk.

# Quick Start Guide
This section shows you how to quickly get set up to factor a number. For a more detailed documentation, see docs/GUIDE.md.

### Set up command machine (e.g., your workstation)
Install Ansible using these [instructions](http://docs.ansible.com/ansible/intro_installation.html#installation)

Set up and configure AWS CLI using these [instructions](http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-set-up.html)

Install GNU Parallel using these [instructions](http://www.gnu.org/software/parallel/)

### Go to the ec2-cluster directory

    >$ cd ec2-cluster

### Configure AWS (Amazon Web Services) environment
The following script will create a new AWS VPC (Virtual Private Cloud) configured for FaaS. 

    >$ ./configure-aws.py

### Edit the values in vars/custom.yml using your favorite text editor (you should be able to figure out what to change)

    >$ vim vars/custom.yml

### Build a base AMI (optional, we have already done this step for you)
We provide a public AMI (ami-e1d39c84) that is the result of running the following script. However, if you wish to create your own base AMI you can run it yourself. It takes 20-30 minutes.

    >$ ./ec2-cluster/build-base.sh

### Run test factorization (optional, but highly recommended)
To check that your AWS environment is correctly configured, and that there are no issues with the setup, we recommend that you run the following test factorization. This following script will build a custom AMI for the test factorization, launch a cluster of four m4.large nodes, and factor a 100-digit number. The entire process should cost less than a dollar in EC2 credit, but will hopefully help you debug any issues with cluster setup. We recommend that you run the commands in the script one by one.

    >$ ./test-factor.sh

### Run factorization
Run the following script to build a new AMI, launch a cluster, and factor the 512-bit key that you specified in vars/custom.yml. We recommend that you run the commands in the script one by one. After you've already built a custom AMI, it is no longer necessary to run the first few commands that build a new AMI.

    >$ ./factor.sh

### Collect results
If all goes to plan, you will recieve an email with the results of the factorization. However, this is research code, and things do not always go as planned :).  

By default, the master node will be stopped (not terminated) after the factorization has completed so that the relevant log files will not be deleted. This will continue to use resources on your Amazon account until the node is terminated. 

The relevant log files are in the following locations on the master node:
    
    /home/ubuntu/server.stderr          # the supervisor output file. You can watch the factorization live with 'tail -f server.stderr'.
    /workdir/<job_name>/<job_name>.log  # the faas log file
    /workdir/<job_name>/<job_name>.cmd  # the commands executed by the factoring script
    /var/log/slurm/slurmctld.log        # the Slurm controller daemon log file

Good luck! 
