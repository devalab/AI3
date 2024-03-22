# Protein Ligand Complex (PLC) Data Generation Workflow on AWS

This repository contains instructions and scripts that are used for creating an AWS Parallel cluster to generate simulation results for 20K PLCs with 5 sample runs for each PLC.

To simplify and automate creation of the compute environment we use a [do-framework](https://bit.ly/do-framework) project called [aws-do-pcluster](https://bit.ly/aws-do-pcluster).

## Build and run the `aws-do-pcluster` container

Clone the [aws-do-pcluster](https://bit.ly/aws-do-pcluster) project. We use this project to build a container which has all required utilities for creation of an Amazon Machine Image and provisioning of the parallel cluster infrastructure.

```bash
git clone https://github.com/aws-samples/aws-do-pcluster
cd aws-do-pcluster
./build.sh
./run.sh
./exec.sh
```
You should now be in a shell within the `aws-do-pcluster` container. To configure access to your AWS account, execute:

```bash
./aws-config.sh
```

## AMI creation

The container comes with a [packer](https://packer.io) recipe that builds an AMI with [GROMACS](https://www.gromacs.org/), [PLUMED](https://plumed.org), and other required modules. 

To create the AMI, simply use the script provided in the container.

```bash
./ami-create.sh
```

An AMI containing all necessary tools and libraries will be created in the configured account and region. 

## VPC creation

The cluster infrastructure requires a VPC. You can list your existing VPCs using the `./vpc-list.sh` script. If you already have a VPC that you'd like to use, list the subnets in your VPC by running `./vps-subnets.sh <vpc_id>`. You will need to select subnets for your head node and your compute nodes. The head node typically resides in a public subnet, while compute nodes use private subnets. 

If you do not have an existing VPC, you can use the [aws-do-cli](https://bit.ly/aws-do-cli) project or other means to create one. 

## Keypair creation

You will use an SSH key to connect to the head node of the cluster, use `./keypair-list.sh` to show a list of existing keypairs in your account. If you do not have a keypair that you would like to use with the head node, to create one execute:

```bash
./keypair-create.sh <keypair_name>
```

where `kaypair_name` is the name of the keypair that you wish to create.

## Cluster setup

### Create Parallel Cluster configuration

From the `aws-do-pcluster` shell, execute:

```bash
./pcluster-config.sh
```

This will open the `pcluster.yaml` file in a `vi` editor.

Copy the content of [pcluster-config.yaml-template](./pcluster-config.yaml-template) into the editor and replace all environment variable placeholders ( $AWS_REGION, $AMI_ID, $KEY_NAME, $SUBNET_ID_HEADNODE, $SUBNET_ID_ONDEMAND, $SUBNET_ID_SPOT_*) with your desired values, using the information from above. Save your changes.

### Create cluster

Once the cluster configuration is created, you are ready to provision the cluster.

```bash
./pcluster-create.sh
```

This creates a cluster named `pcluster` from the configuration file named `pcluster.yaml`. You can specify a different cluster name or configuration file by passing arguments to the script. Example: `./pcluster-create.sh <cluster_name> <config_name>` 

The cluster creation process takes approx. 20 minutes.

### Connect to the cluster

Once the cluster is created, you can connect to the head node by using the `./pcluster-connect.sh` script.

```bash
./pcluster-connect.sh [cluster_name] [kay_path]
```

If executed without arguments, the script will use `pcluster` as cluster name and `/wd/ssh/pcluster.pem` as the key_path.

## Array job submissions

The cluster has a shared FSx for Lustre volume mounted in path `/fsx`. From the head node, clone this repository in the shared path `/fsx`, so the scripts can be accessible on both the head node and compute nodes. You will be executing all further commands from the head node shell.

1) `./activate_environment.sh` - activate virtual environment for running and processing PLC jobs on the cluster. 

2)	`./submit_dataset_array.sh` - script to setup environment variables dynamically for job submission with configurations from the .env-template file. 

Usage: `./submit_dataset_array.sh <start_index> <end_index> <batchNumber> <partition> <coreNumber> <skipFailedPLC>`

E.g., 
```
./submit_dataset_array.sh 12000 13999 VII large 24| tee -a /tmp/submit_logs/large7
```

3)	`./submit_dataset_array.py` - Python script to run a batch of jobs – both the finished or submitted jobs will be skipped

4)	`./run_scripts_array.sh` - submit an array job for individual PLC

5)	submit_plas_array.sh: a template file used to create a script for array job submission at runtime

6)	run_simulations.py: Python script to run individual simulation 

## Result summary
7)	count_results.sh: Summarize the total number of PLC jobs finished successfully across multiple batches set in the count_locations.conf file

## Sync the result from a cluster to Amazon S3
8)	sync_results-dynamic.sh: Sync the final result to Amazon S3 for a batch to make sure all successful PLC results uploaded

Usage: `sync_results-dynamic.sh <batch number>`


As a result of this workflow, you will have the generated PLC data stored in your S3 bucket.

## Cleanup

The Parallel Cluster infrastructure can be removed by simply executing the delete script from the `aws-do-pcluster` container shell:

```bash
./pcluser-delete.sh
```

The generated data that is stored in your S3 bucket will not be deleted.


## Security

See [CONTRIBUTING](CONTRIBUTING.md#security-issue-notifications) for more information.

## License

This library is licensed under the MIT-0 License. See the LICENSE file.
