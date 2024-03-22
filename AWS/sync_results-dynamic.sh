#!/bin/bash

#source .env
S3_OUTPUT_PATH=s3://plcrun/output-20K-$1
S3_LOG_PATH=s3://plcrun/logs-20K-$1
OUTPUT_PATH=/fsx/plcrun/output-20K-$1
LOG_PATH=/fsx/plcrun/logs-20K-$1
NUM_RUNS=5

pushd $OUTPUT_PATH

ls -d */ | sed 's/\///g' > dir.txt

cnt=0
for d in $(cat dir.txt) ; do
	echo "Processing output directory ${d} ..."
	finished_num=$(find ${d} -name "FINAL_RESULTS*" | wc -l)
	if [ "$NUM_RUNS" == "$finished_num" ]; then
		# Sync output
		echo "All runs for PLC ${d} have finished. Synching ..."
		tar czvf ${d}.tar.gz ${d}
		aws s3 cp ${d}.tar.gz ${S3_OUTPUT_PATH}/${d}.tar.gz
		rm ${d}.tar.gz
		# Sync logs
		pushd $LOG_PATH
		tar czvf ${d}-logs.tar.gz ${d}
		aws s3 cp ${d}-logs.tar.gz ${S3_LOG_PATH}/${d}-logs.tar.gz
		rm ${d}-logs.tar.gz
		popd
	else
		echo "$finished_num of $NUM_RUNS runs for PLC ${d} have finished. Skipping sync ..."
	fi

	cnt=$((cnt+1))
done

popd

echo ""
echo "Processed $cnt directories"
echo ""

