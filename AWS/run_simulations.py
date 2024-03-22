import argparse
import logging
import os
import subprocess
import sys
import time

from dotenv import load_dotenv

# Paths
load_dotenv()

PARAM_PATH = os.getenv("PARAM_PATH")
INPUT_PATH = os.getenv("INPUT_PATH")
OUTPUT_PATH = os.getenv("OUTPUT_PATH")
LOG_PATH = os.getenv("LOG_PATH")

BASH = "/bin/bash"
CPT_TIME = os.getenv("CPT_TIME")  # Checkpointing time
NUM_PROC = os.getenv("NUM_PROC")  # Number of cores to be used for the run

# NVT simulation: Heating of system
def run_nvt() -> None:
    command = f"""
        gmx grompp -f {PARAM_PATH}/nvt.mdp -o nvt.tpr -c {INPUT_PATH}/em.gro -r {INPUT_PATH}/em.gro -p {INPUT_PATH}/complex_solvated_prot.top -n {INPUT_PATH}/indx.ndx -maxwarn 1
    """
    subprocess.run(command, shell=True, check=True, executable=BASH)

    try:
        command = f"gmx mdrun -nt {NUM_PROC} -deffnm nvt -cpt {CPT_TIME} -cpi nvt.cpt -v"
        subprocess.run(command, shell=True, check=True, executable=BASH)
    except:
        command = f"gmx mdrun -nt {NUM_PROC} -deffnm nvt -cpt {CPT_TIME} -cpi nvt_prev.cpt -v"
        subprocess.run(command, shell=True, check=True, executable=BASH)


# Equilibration and Production Run
def run_prod() -> None:
    command = f"""
        gmx grompp -f {PARAM_PATH}/prod.mdp -o prod.tpr -c nvt.gro -p {INPUT_PATH}/complex_solvated_prot.top -maxwarn 1
    """
    subprocess.run(command, shell=True, check=True, executable=BASH)

    try:
        command = f"gmx mdrun -nt {NUM_PROC} -deffnm prod -cpt {CPT_TIME} -cpi prod.cpt -v"
        subprocess.run(command, shell=True, check=True, executable=BASH)
    except:
        command = f"gmx mdrun -nt {NUM_PROC} -deffnm prod -cpt {CPT_TIME} -cpi prod_prev.cpt -v"
        subprocess.run(command, shell=True, check=True, executable=BASH)


# Post-production Analysis
def mmpbsa() -> None:
    group = []

    with open(f'{INPUT_PATH}/indx.ndx','r') as f:
        for line in f:
            if "[" in line:
                group.append(line.rstrip())

    grp = group.index('[ LIG ]')

    command = f"""
        echo {grp} 0 | gmx trjconv -s prod.tpr -f prod.xtc -o prod_cent.xtc -pbc mol -center

        echo 1 1 | gmx rms -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsd_1.xvg -mw

        echo 3 3 | gmx rms -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsd_3.xvg -mw

        echo 4 4 | gmx rms -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsd_4.xvg -mw

        echo 5 5 | gmx rms -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsd_5.xvg -mw

        echo 13 13 | gmx rms -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsd_13.xvg -mw

        echo 14 14 | gmx rms -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsd_14.xvg -mw

        echo 3 | gmx rmsf -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsf_3.xvg -od rmsfd_3.xvg

        echo 4 | gmx rmsf -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsf_4.xvg -od rmsfd_4.xvg

        echo 5 | gmx rmsf -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsf_5.xvg -od rmsfd_5.xvg

        echo 13 | gmx rmsf -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsf_13.xvg -od rmsfd_13.xvg

        echo 14 | gmx rmsf -s prod.tpr -f prod_cent.xtc -n {INPUT_PATH}/indx.ndx -o rmsf_14.xvg -od rmsfd_14.xvg

        if [ -f "{INPUT_PATH}/cofactor_prot_0.pdb" ]; then mpirun -np {NUM_PROC} gmx_MMPBSA -nogui -O -i {PARAM_PATH}/mmpbsa.in -cs prod.tpr -ci {INPUT_PATH}/indx.ndx -eo Interaction_energy.csv -cg 1 14 -lm {INPUT_PATH}/ligand_clean_h_prot.mol2 -ct prod_cent.xtc > progress.log
        else
        mpirun -np {NUM_PROC} gmx_MMPBSA -nogui -O -i {PARAM_PATH}/mmpbsa.in -cs prod.tpr -ci {INPUT_PATH}/indx.ndx -eo Interaction_energy.csv -cg 1 13 -lm {INPUT_PATH}/ligand_clean_h_prot.mol2 -ct prod_cent.xtc > progress.log
        fi
    """
    subprocess.run(command, shell=True, check=True, executable=BASH)


def is_process_done(run_type: str) -> bool:
    """Check if the process is done or not

    Args:
        run_type (str): nvt|prod

    Returns:
        bool: True if the process is done, False otherwise
    """

    if not os.path.exists(f"{run_type}.log"):
        return False
    output = subprocess.check_output(
        f"tail -n 3 {run_type}.log", shell=True, encoding="utf-8"
    ).split("\n")
    for out in output:
        if "Finished mdrun" in out:
            return True
    return False


def clean_dir() -> None:
    """Clean the directory of all the unnecessary files"""

    command = "rm --force \#* _* || true"
    subprocess.run(command, shell=True, check=True, executable=BASH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scripts")
    parser.add_argument("index", type=str, nargs="?")
    parser.add_argument("pdbid", type=str, nargs="?")
    parser.add_argument("run_idx", type=int, nargs="?")
    parser.add_argument("-l", "--log", action="store_true")

    log = parser.parse_args().log

    if parser.parse_args().index and parser.parse_args().pdbid:
        args = parser.parse_args()
        index = args.index
        pdbid = args.pdbid
        run_idx = args.run_idx
    else:
        print("Please enter the index, pdbid and run_idx")
        sys.exit(1)

    INPUT_PATH = f"{INPUT_PATH}/{index}.{pdbid}"
    OUTPUT_PATH = f"{OUTPUT_PATH}/{index}.{pdbid}"
    LOG_PATH = f"{LOG_PATH}/{index}.{pdbid}/{index}_{pdbid}_{run_idx}.log"

    if log:
        logging.basicConfig(
            filename=f"{LOG_PATH}",
            filemode="a",
	    level=logging.INFO,
            format="[%(asctime)s]: [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        logging.Formatter.converter = time.gmtime

    logging.info("Running job %s.%s", index, pdbid)

    run_path = f"{OUTPUT_PATH}/{index}.{pdbid}_{run_idx}"

    if not os.path.exists(run_path):
        logging.info("Creating directory %s", run_path)
        os.makedirs(run_path)
    os.chdir(run_path)

    if not is_process_done("nvt"):
        logging.info("Running NVT simulation")
        try:
            run_nvt()
            logging.info("NVT simulation done")
        except subprocess.CalledProcessError as err:
            logging.error("NVT simulation failed")
            logging.error(err)
            sys.exit(1)
    else:
        logging.info("NVT simulation already done")

    if not is_process_done("prod"):
        logging.info("Performing Equilibration and Production Run")
        try:
            run_prod()
            logging.info("Equilibration and Production Run done")
        except subprocess.CalledProcessError as err:
            logging.error("Equilibration and Production Run failed")
            logging.error(err)
            sys.exit(1)
    else:
        logging.info("Equilibration and Production Run already done")

    if not os.path.exists(f"{run_path}/FINAL_RESULTS_MMPBSA.dat"):
        logging.info("Running MMPBSA")
        try:
            mmpbsa()
            logging.info("MMPBSA done")
        except subprocess.CalledProcessError as err:
            logging.error("MMPBSA failed")
            logging.error(err)
            sys.exit(1)
    else:
        logging.info("MMPBSA already done")

    try:
        clean_dir()
    except subprocess.CalledProcessError as err:
        logging.error(err)
        sys.exit(1)

    logging.info("All processes done")
