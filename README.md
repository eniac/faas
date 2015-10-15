# FaaS
The purpose of the FaaS (Factoring as a Service) project is to demonstrate that 512-bit integers can be factored in only a few hours, for less than $100 of compute time in a public cloud environment.  This illustrates the amazing progress in computing power over time, and the risk of continued use of 512-bit RSA keys.

Please do not use these scripts to attack any systems that you do not own.

Our scripts launch a compute cluster on Amazon EC2 and run the [CADO-NFS](http://cado-nfs.gforge.inria.fr/) and [Msieve](http://sourceforge.net/projects/msieve/) implementations of the number field sieve factoring algorithm, with some improvements to parallelization. For more information about the project, see [the project webpage](https://seclab.upenn.edu/projects/faas/) and [our paper](https://seclab.upenn.edu/projects/faas/faas.pdf).

The tool was built using the AWS region 'us-east' and currently has many hard-coded settings for this region.

# Quick Start Guide
This section shows you how to quickly get set up to factor. 

### Set up command machine (e.g., your workstation)
Set up and configure AWS CLI using these [instructions](http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-set-up.html).
Make sure that your ~/.aws/config looks like the following:

    [default]
    region = <EC2-region>
    output = json

and your ~/.aws/credentials looks like

    [default]
    aws_access_key_id = <key_id>
    aws_secret_access_key = <access_key>

and your ~/.boto config looks like

    [Credentials]
    aws_access_key_id = <key_id>
    aws_secret_access_key = <access_key>

Install Ansible using these [instructions](http://docs.ansible.com/ansible/intro_installation.html#installation) and configure Boto (python interface to AWS) using these [instructions](http://docs.ansible.com/ansible/intro_dynamic_inventory.html#example-aws-ec2-external-inventory-script).

Install GNU Parallel using these [instructions](http://www.gnu.org/software/parallel/).

### Go to the ec2 directory (some of our scripts use relative paths)

    >$ cd ec2

### Download Ansible EC2 dynamic inventory scripts
The scripts are available [here](http://docs.ansible.com/ansible/intro_dynamic_inventory.html#example-aws-ec2-external-inventory-script), and you can download them with the following commands:

    >$ wget https://raw.githubusercontent.com/ansible/ansible/devel/contrib/inventory/ec2.py
    >$ chmod +x ec2.py
    >$ wget https://raw.githubusercontent.com/ansible/ansible/devel/contrib/inventory/ec2.ini

Set the following values in ec2.ini.

    regions = <EC2-region>  # you can leave this as 'all', but the script runs more quickly with a single region
    cache_max_age = 0       # we want to always see the most up-to-date info about running instances, so do not cache
    rds = False             # if your AWS user does not have rds permissions
    elasticache = False     # if your AWS user does not have elasticache permissions

### Configure AWS (Amazon Web Services) environment
The following script will create a new AWS VPC (Virtual Private Cloud) configured for FaaS. 

    >$ ./configure-aws.py

### Set custom values for your setup in vars/custom.yml using your favorite text editor.

    >$ vim vars/custom.yml

### Build a base AMI
We provide a public AMI for the region us-east-1 (ami-19642b7c), which is the result of running the following script. To use our public AMI, replace the value of the 'base\_image' variable in 'vars/ec2.py' with 'ami-19642b7c' and do not run the script. The script can take up to an hour to run.

    >$ ./build-base.sh

### Optional: Run test factorization
To check that your AWS environment is correctly configured, we recommend that you run a small test factorization. Our test script will build a custom AMI for the test factorization, launch a cluster of four m4.large nodes, and factor a 100-digit number. The entire process should cost less than a dollar in EC2 credit, but will hopefully help you debug any issues with cluster setup. We recommend that you run the commands in the script one by one.

    >$ ./test-factor.sh

### Run factorization
The following script will build a new AMI, launch a cluster, and factor the 512-bit integer that you specified in vars/custom.yml. We recommend that you run the commands in the script one by one. After you've already built a custom AMI, it is no longer necessary to build a new AMI each time, so you may wish to comment out the first few lines.

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
