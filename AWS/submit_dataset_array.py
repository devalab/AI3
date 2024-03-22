import sys
import os
from glob import glob

from dotenv import load_dotenv

start = int(sys.argv[-2])
end = int(sys.argv[-1])

# Counters
cnt_job_new=0
cnt_job_existing=0
cnt_job_completed=0
cnt_plc_total=end-start+1
cnt_plc_inrange=0
cnt_plc_outrange=0

# Paths
load_dotenv()

INPUT_PATH = os.getenv("INPUT_PATH")
STORE_PATH = os.getenv("STORE_PATH")
NUM_RUNS = int(os.getenv("NUM_RUNS"))
SCRIPT_PATH = os.getenv("SCRIPT_PATH")
batchNumber = os.getenv("batchNumber")

job_list=""
with open(f'{SCRIPT_PATH}/job_queue_array-{batchNumber}.txt', 'r') as job_list_file:
        job_lines = job_list_file.readlines()
        job_list=' '.join(job_lines)

x = glob(f"{INPUT_PATH}/*", recursive = False)
x.sort()

for i in range(start, end+1):
    folder = x[i].split("/")[-1]
    index = folder.split(".")[0]
    pdbid = folder.split(".")[1]

    idxlist=[]
    for r in range(0,NUM_RUNS):
        finished = os.path.exists(f'{STORE_PATH}/{index}.{pdbid}/{index}.{pdbid}_{r}/FINAL_RESULTS_MMPBSA.dat')
        if (finished):
            cnt_job_completed = cnt_job_completed + 1
            print(f'Job {index}.{pdbid}_{r} is already completed. Skipping ...')
            continue
        else:
            # Check if job is already submitted
            job_name=f'{index}_{pdbid}_{r}'
            if job_name in job_list:
                cnt_job_existing = cnt_job_existing + 1
                print(f'Job {job_name} already exists. Skipping submission ...')
                sys.stdout.flush()
            else:
                cnt_job_new = cnt_job_new + 1
                print(f'Job {job_name} has never been run. Adding to the list ...')
                sys.stdout.flush()
                idxlist.append(str(r))
    if (len(idxlist) != 0):
        idxStr = ','.join(idxlist)
        os.system(f'bash run_scripts_array.sh {index} {pdbid} {idxStr}')

print('\nJob submission summary:')
print(f'Start index: {start}')
print(f'End index: {end}')
print(f'Total PLCs: {cnt_plc_total}')
print(f'Jobs already completed: {cnt_job_completed}')
print(f'Jobs in progress: {cnt_job_existing}')
print(f'Jobs new: {cnt_job_new}')

