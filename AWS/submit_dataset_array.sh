#!/bin/bash

if [ "$5" == "" ]; then                                                                           
	echo ""                                                                                    
	echo "Usage:"                                                                              
	echo "   $0 <start_index> <end_index> <batchNumber> <partition> <coreNumber> <skipFailedPLC>"                 
	echo ""                                                                                     
else 
	source activate gmx_MMPBSA

	export batchNumber=$3
	export instanceSize=$4
	export coreNumber=$5
        skipFailedPLC=$6
fi

export $(grep -v '^#' .env-template | sed -e "s/batchNumber/$batchNumber/" -e "s/instanceSize/$instanceSize/" -e "s/coreNumber/$coreNumber/" | xargs -d '\n') 

date > $SCRIPT_PATH/job_queue_array-$batchNumber.txt
squeue -O name,jobarrayid |awk -F " +|\t|-" '{name=$1; if (NF >3)  { gsub(".*_\\[", "", $2); gsub ("\\]", "", $3); for(i=$2;i<=$3;i++) print $1 "_" i} else if ($2 ~ "_") {gsub(".*_","",$2); print $1 "_" $2} else print $1 }' >> $SCRIPT_PATH/job_queue_array-$batchNumber.txt 

# Skip completely failed PLCs
if [[ $skipFailedPLC == "skip" ]]; then
  pushd /fsx/plcrun/output-20K-$batchNumber/
  for d in `\ls `; do number=`find $d -name FINAL_RESULTS_MMPBSA.dat| wc -l`; if [[ $number == 0 ]]; then echo ${d/./_}_{0..4}; fi ; done >> $SCRIPT_PATH/job_queue_array-$batchNumber.txt
  popd
fi

python submit_dataset_array.py $1 $2

