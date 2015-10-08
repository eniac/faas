# FaaS
The purpose of the FaaS (Factoring as a Service) project is to create a tool that allows anyone to factor a 512-bit RSA key with minimal time and cost. To do this, we have developed a framework for launching an compute cluster on Amazon EC2 and running a parallelized implementation of the Number Field Sieve (NFS) algorithm to factor keys.

The tool was built using the AWS region 'us-east'. Run in any other region at your own risk.

# Quick Start
This section shows you how to quickly get set up to factor a number. For a more detailed documentation, see docs/GUIDE.md.

1.  Set up command machine (e.g., your workstation)
Install Ansible using these [instructions](http://docs.ansible.com/ansible/intro_installation.html#installation)

Set up and configure AWS CLI using these [instructions](http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-set-up.html)

2. Configure AWS (Amazon Web Services) environment
The following script will create a new AWS VPC (Virtual Private Cloud) configured for FaaS. 
    >$ ./configure-aws.py

3. Create custom AMI (Amazon Machine Image)
In order to facilitate a quick setup process, we provide a public AMI with default configurations. This AMI was build using the `build-base.sh` script. However, you will need a custom AMI with newly-generated SSH and Munge keys, along with any other custom configurations that you desire. Run the following script to build a custom AMI.

    >$ ./build-custom.sh

4. Run test factorization (optional, but highly recommended)
To check that your AWS environment is correctly configured, and that there are no issues with the setup, we recommend that you run the following test factorization. This will launch a cluster of four TODO nodes to factor a 100-digit number. The entire process should take less than 20 minutes and should cost only a few dollars in EC2 credit.

    >$ ./test-factor.sh

5. Run factorization
Once the previous steps have been completed, it should be possible to start a new factorization job only varying the key to be factored. Edit ec2\_variables.yml and change params.N to the number (in decimal) that you wish to factor.

    >$ ./factor.sh

6. Collect results
If you run the factorization using supervisor, results should be emailed to you. If you do not see an email, you can SSH into the master node to check the log files.

By default, the master node will be stopped (not terminated) after the factorization has completed. This will continue to use resources on your Amazon account until the node is terminated. 
