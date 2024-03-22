#!/bin/bash

#. .env-test

index=$1
pdbid=$2
runlst=$3

num_cores=$NUM_CORES_PER_NODE #Total available maximum number of cores at the node.
account=d4
partition=$PARTITION_NAME
num_runs=$NUM_RUNS #Number of independent runs.

# paths
#export $(grep -v '^#' .env-test | xargs -d '\n')

input_dir="$INPUT_PATH"
script_dir="$SCRIPT_PATH"
store_dir="$STORE_PATH/${index}.${pdbid}/${index}.${pdbid}"
log_path="$LOG_PATH/${index}.${pdbid}"

mkdir -p $log_path

for (( run_idx=0; run_idx<$num_runs; run_idx++ ))
do
    mkdir -p "${store_dir}_${run_idx}" 
done   

submit_script="${script_dir}/${index}_${pdbid}.sh"
log_path_curr="${log_path}/${index}_${pdbid}_\%a.out"
cp ${script_dir}/submit_plas_array.sh $submit_script
sed -i "s/pdbid/${pdbid}/g" $submit_script
sed -i "s/index/${index}/g" $submit_script
sed -i "s/account/${account}/g" $submit_script
sed -i "s/partition/${partition}/g" $submit_script
sed -i "s/job_name/${index}_${pdbid}/g" $submit_script
sed -i "s/num_cores/${num_cores}/g" $submit_script
sed -i "s/run_idx/${run_idx}/g" $submit_script
sed -i "s#log_path#${log_path_curr}#g" $submit_script
sbatch --array=$runlst ${submit_script}
rm ${submit_script}

