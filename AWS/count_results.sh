#!/bin/bash
NUM_RUNS=5
total_finished=0
total_partial=0
total_finished_jobs=0
total_dir=0
locations=0

for o in $(cat ./count_locations.conf); do

	echo ""
	echo "Processing batch location: $o ..."
	locations=$((locations+1))
	pushd $o

	ls -d */ | sed 's/\///g' > dir.txt

	finished_cnt=0
	partial_cnt=0
	finished_jobs_cnt=0
	cnt=0
	for d in $(cat dir.txt) ; do
		#echo "Processing output directory ${d} ..."
		finished_num=$(find ${d} -name "FINAL_RESULTS_MMPBSA.dat" | wc -l)
		finished_jobs_cnt=$((finished_jobs_cnt+finished_num))
		if [ "$NUM_RUNS" == "$finished_num" ]; then
			#echo "All runs for PLC ${d} have finished. Adding to finished count ..."
			finished_cnt=$((finished_cnt+1))
		else
			if [ $finished_num -gt 0 ]; then
				partial_cnt=$((partial_cnt+1))
			fi
		fi

		cnt=$((cnt+1))
	done

	echo ""
	echo "Batch location: $o"
        echo "Submitted PLCs: $cnt"
	echo "Submitted Jobs: $((cnt*NUM_RUNS))"
	echo "Finished Jobs: $finished_jobs_cnt"
        echo "Finished PLCs: $finished_cnt"
	echo "Partial PLCs: $partial_cnt"
	echo "Total PLCs: $((finished_cnt+partial_cnt))"
	total_dir=$((total_dir+cnt))
	total_finished=$((total_finished+finished_cnt))
	total_partial=$((total_partial+partial_cnt))
	total_finished_jobs=$((total_finished_jobs+finished_jobs_cnt))

	popd

done

echo ""
echo "Total batches: $locations"
echo "Total Submitted PLCs: $total_dir"
echo "Total Submitted Jobs: $((total_dir*NUM_RUNS))"
echo "Total Finished Jobs: $((total_finished_jobs))"
echo "Total Finished PLCs: $total_finished"
echo "Total Partial PLCs: $total_partial"
echo "Grand Total PLCs: $((total_finished+total_partial))"
echo ""

