import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import tarfile

from parmed.exceptions import FormatNotFound
from prody import parsePDB, writePDB
from termcolor import colored

from scripts.utils import (
    amber_to_gmx,
    download_protoss_pdb,
    file_map,
    get_atom_count,
    get_charge,
    itp_reindexing,
    open_ligand_in_browser,
    parse_pdb,
)

BASH = "/bin/bash"
ZSH = "/bin/zsh"

# Paths
absolute_path = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.join(absolute_path, "data")
input_path = os.path.join(data_path, "input")
work_path = os.path.join(data_path, "work")
script_path = os.path.join(absolute_path, "scripts")
download_path = os.path.join(data_path, "download")
output_path = os.path.join(data_path, "output")
conda_env_path = os.path.dirname(os.path.dirname(sys.executable))
ff_location = os.path.join(conda_env_path, "dat", "leap", "cmd")

if not os.path.exists(ff_location):
    raise FileNotFoundError(
        f"Could not find the FF files at {ff_location}. Please check if the conda environment is activated."
    )

index, pdbid = "", ""
cofactors = []
cofactor_charges = []
ligand_charge = []


def prepare_input():
    print("Preparing input files...")
    tar_file = f"{index}.{pdbid}.tar"
    src = os.path.join(input_path, tar_file)
    dst = os.path.join(work_path, tar_file)
    shutil.copy(src, dst)

    with tarfile.open(dst) as tar:
        tar.extractall(path=work_path)
    os.remove(dst)

    complex_solvated_path = os.path.join(
        work_path, f"{index}.{pdbid}", "complex_solvated.pdb"
    )
    if not os.path.isfile(complex_solvated_path):
        raise FileNotFoundError

    print()


def clean_input():
    pdbfile = os.path.join(download_path, f"{pdbid}.pdb")
    if not os.path.exists(pdbfile):
        download_protoss_pdb(pdbid, download_path)

    cof = glob.glob("cofactor*pdb")
    print()
    print("--- Processing cofactors ---")
    print(
        colored(
            f"#Info: There are {len(cof)} cofactors in this complex.",
            "green",
            attrs=["bold"],
        )
    )
    for i, cofactor in enumerate(cof):
        print(f"Processing cofactor {i + 1}...")
        cf, chain_id, res_idx = parse_pdb(cofactor)
        cofactors.append(cf)
        print(
            colored(
                f"# Info: Chain ID - {chain_id}, Residue Index - {res_idx}",
                "green",
                attrs=["bold"],
            )
        )
        parsed_cofactor = parsePDB(pdbfile)
        output = f"cofactor_prot_{i}.pdb"
        cofact = parsed_cofactor.select(
            f"resname {cf} and chain {chain_id} and resid {res_idx}"
        )
        writePDB(f"cofactor_prot_{i}.pdb", cofact)
        os.system(f"""sed -i "s/\'/Z/g" cofactor_prot_{i}.pdb""")
        try:
            predicted_charge = get_charge(
                f"{work_path}/cofactor_prot_{i}.pdb",
                f"{script_path}/get_charge.sh",
            )
            print(
                colored(
                    f"# Info: Chimera predicted charge: {predicted_charge}",
                    "green",
                    attrs=["bold"],
                )
            )
        except (subprocess.CalledProcessError, UnboundLocalError) as err:
            print(colored(err, "red", attrs=["bold"]))

        if not is_automated:
            open_ligand_in_browser(cofactors[i])
            subprocess.run(
                f"chimera cofactor_prot_{i}.pdb",
                shell=True,
                check=True,
                executable=BASH,
            )

            cofactor_charges.append(
                int(input(f"Enter the charge of cofactor {i + 1}: "))
            )
        else:
            cofactor_charges.append(predicted_charge)

    print()
    print("--- Processing ligand ---")
    if os.path.exists("ligand.pdb"):
        _, chain_id, res_idx = parse_pdb("ligand.pdb")
        print(
            colored(
                f"# Info: Chain ID - {chain_id}, Residue Index - {res_idx}",
                "green",
                attrs=["bold"],
            )
        )
        mol_type = "lig"
        file_map("ligand.pdb", "ligand_prot.pdb", mol_type, pdbfile)
        output = "ligand_prot.pdb"
        ligand_id, *_ = parse_pdb(output)

        with open(output, "r", encoding="utf-8") as f:
            for _ in range(5):
                print(f.readline().strip())
            print(colored(f"# Info: Ligand ID is {ligand_id}", "green", attrs=["bold"]))

        try:
            predicted_charge = get_charge(
                f"{work_path}/ligand_prot.pdb",
                f"{script_path}/get_charge.sh",
            )
            print(
                colored(
                    f"# Info: Chimera predicted charge: {predicted_charge}",
                    "green",
                    attrs=["bold"],
                )
            )
        except (subprocess.CalledProcessError, UnboundLocalError) as err:
            print(colored(err, "red", attrs=["bold"]))

        if not is_automated:
            open_ligand_in_browser(ligand_id)
            subprocess.run(
                "chimera ligand_prot.pdb", shell=True, check=True, executable=BASH
            )

            ligand_charge.append(int(input("Enter the charge of the ligand: ")))
        else:
            ligand_charge.append(predicted_charge)

    elif os.path.exists("aminoacid.pdb"):
        mol_type = "aminoacid"
        file_map("aminoacid.pdb", "aminoacid_prot.pdb", mol_type, pdbfile)
    elif os.path.exists("peptide.pdb"):
        mol_type = "peptide"
        file_map("peptide.pdb", "peptide_prot.pdb", mol_type, pdbfile)

    return mol_type


def generate_params():
    command = f"sh {script_path}/second_pargen.sh"
    num_cofactors = len(cofactors)
    shell_input = f"{index}\n{pdbid}\n{num_cofactors}\n"
    shell_input += f"{ligand_charge[0]}\n"
    if num_cofactors != 0:
        for i in range(num_cofactors):
            shell_input += f"{cofactor_charges[i]}\n"

    subprocess.run(
        command,
        input=shell_input,
        encoding="utf-8",
        shell=True,
        check=True,
        executable=BASH,
    )


def combine_all():
    command = f"sh {script_path}/third_combineall.sh {index} {pdbid}"
    subprocess.run(
        command,
        encoding="utf-8",
        shell=True,
        check=True,
        executable=BASH,
    )


def convert_to_gromacs():
    print()
    print("--- Converting to gromacs... ---")
    command = f"cp {script_path}/*mdp ."
    subprocess.run(command, shell=True, check=True, executable=BASH)
    return amber_to_gmx(tuple(cofactors))


def count_atoms():
    ligand_count, cofactor_count = [], []
    ligand_protoss = get_atom_count("ligand_prot.pdb")
    ligand_ada = get_atom_count("ligand.pdb")
    ligand_count.append(f"{ligand_protoss}-{ligand_ada}")
    for i in range(len(cofactors)):
        cofactor_protoss = get_atom_count(f"cofactor_prot_{i}.pdb")
        cofactor_ada = get_atom_count(f"cofactor_{i}.pdb")
        cofactor_count.append(f"{cofactor_protoss}-{cofactor_ada}")
    return ligand_count, cofactor_count


def print_data():
    global num_atoms
    num_atoms = int(
        subprocess.check_output("sed -n 2p complex_solvated_prot.gro", shell=True)
    )

    print(colored("Data for the excel sheet", "green", attrs=["bold"]))
    print(colored(f"index: {index}", "green", attrs=["bold"]))
    print(colored(f"pdbid: {pdbid}", "green", attrs=["bold"]))
    print(colored(f"Ligand Charge: {ligand_charge}", "green", attrs=["bold"]))
    print(colored(f"Cofactor Charge: {cofactor_charges}", "green", attrs=["bold"]))
    print(colored(f"Number of atoms: {num_atoms}", "green", attrs=["bold"]))


def upload():
    print_data()
    command = "rm "
    if os.path.exists("index.ndx"):
        command += "index.ndx "
    for file in glob.glob("#*"):
        command += f"'{file}' "
    if os.path.exists("geckodriver.log"):
        command += "geckodriver.log"
    subprocess.run(command, shell=True, check=True, executable=BASH)
    print()
    print(colored("Do you want to upload the files to ada?", "yellow", attrs=["bold"]))
    print(
        colored(
            "(Perform this step only if you are sure everything is correct)",
            "yellow",
            attrs=["bold"],
        )
    )
    if is_automated:
        inp = "y"
    else:
        inp = input("Enter Y/N: ")
    if inp.strip().lower() == "y":
        os.chdir(absolute_path)
        command = f"tar -czf {index}.{pdbid}.tar {index}.{pdbid}"
        subprocess.run(command, shell=True, check=True, executable=BASH)
        command = f"rm -r {index}.{pdbid}"
        # subprocess.run(command, shell=True, check=True, executable=BASH)
        command = f"""
        sshpass -p "password" scp {index}.{pdbid}.tar d4@ada.iiit.ac.in:{output_path}
        """
        subprocess.run(command, shell=True, check=True, executable=BASH)
        print("Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scripts")
    parser.add_argument("index", type=str, nargs="?")
    parser.add_argument("pdbid", type=str, nargs="?")
    parser.add_argument("-r", "--read_input", action="store_true")
    parser.add_argument("-a", "--automate", action="store_true")
    parser.add_argument("--steps", type=str, default="0-8")
    first_step, last_step = parser.parse_args().steps.split("-")
    steps = tuple(range(int(first_step), int(last_step) + 1))
    is_automated = parser.parse_args().automate
    if parser.parse_args().read_input:
        input_list = f"{input_path}/input_list.txt"
    elif parser.parse_args().index and parser.parse_args().pdbid:
        args = parser.parse_args()
        index = args.index
        pdbid = args.pdbid
        input_list = f"{data_path}/input_list.txt"
        with open(input_list, "w", encoding="utf-8") as f:
            f.write(f"{index} {pdbid}")
    else:
        print("Please enter the index and pdbid, or provide an input file")
        sys.exit()

    out_log = f"{data_path}/output_log.csv"
    already_exists = os.path.exists(out_log)
    with open(input_list, "r", encoding="utf-8") as id_list, open(
        out_log, "a", encoding="utf-8"
    ) as final_out:
        writer = csv.writer(final_out)
        if not already_exists:
            writer.writerow(
                [
                    "index",
                    "pdbid",
                    "Ligand Charge",
                    "Cofactor Charge",
                    "Ligand Count (Protoss-Ada)",
                    "Cofactor Count (Protoss-Ada)",
                    "Atom Count",
                    "Remarks",
                ]
            )
        for line in id_list:
            if not line.strip():
                continue
            index, pdbid = line.strip().split()
            if not index or not pdbid:
                continue

            cofactors = []
            cofactor_charges = []
            ligand_charge = []
            num_atoms = 0

            print("Here are the steps you can perform:")
            print("0. Prepare input files")
            print("1. Clean input files")
            print("2. Generate params")
            print("3. Combine all the files")
            print("4. Convert to gromacs")
            print("5. ITP reindexing")
            print("6. Check in VMD")

            # Zeroth, prepare input files
            try:
                if 0 in steps:
                    prepare_input()

            except subprocess.CalledProcessError as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [
                        index,
                        pdbid,
                        "Error",
                        "Error",
                        "Error",
                        "Error",
                        "Error",
                        "Error (Input file not present)",
                    ]
                )
                continue
            except FileNotFoundError as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [
                        index,
                        pdbid,
                        "Error",
                        "Error",
                        "Error",
                        "Error",
                        "Error",
                        "Error (complex_solvated.pdb not found)",
                    ]
                )
                continue
            except Exception as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [index, pdbid, "Error", "Error", "Error", "Error (Step 0)"]
                )
                continue

            work_path = os.path.join(work_path, f"{index}.{pdbid}")
            os.chdir(work_path)

            print(colored("Step 0 completed successfully", "green", attrs=["bold"]))

            # First, clean the input files
            try:
                if 1 in steps:
                    if not is_automated:
                        input("Press enter to continue...")
                    clean_input()
            except UnboundLocalError as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [
                        index,
                        pdbid,
                        "Error",
                        "Error",
                        "Error",
                        "Error",
                        "Error",
                        "Error (Chimera was unable to predict charge)",
                    ]
                )
                continue
            except Exception as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [index, pdbid, "Error", "Error", "Error", "Error (Step 1)"]
                )
                continue

            print(colored("Step 1 completed successfully", "green", attrs=["bold"]))

            ligand_count, cofactor_count = "", ""
            try:
                ligand_count, cofactor_count = count_atoms()
            except Exception as err:
                print(colored(err, "red", attrs=["bold"]))

            # Second, generate params
            try:
                if 2 in steps:
                    if not is_automated:
                        input("Press enter to continue...")
                    generate_params()
            except Exception as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [
                        index,
                        pdbid,
                        ligand_charge,
                        cofactor_charges,
                        ligand_count,
                        cofactor_count,
                        "Error",
                        "Error (Step 2)",
                    ]
                )
                continue

            print(colored("Step 2 completed successfully", "green", attrs=["bold"]))

            # Third, combine all the files
            try:
                if 3 in steps:
                    if not is_automated:
                        input("Press enter to continue...")
                    combine_all()
            except Exception as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [
                        index,
                        pdbid,
                        ligand_charge,
                        cofactor_charges,
                        ligand_count,
                        cofactor_count,
                        "Error",
                        "Error (Step 3)",
                    ]
                )
                continue

            print(colored("Step 3 completed successfully", "green", attrs=["bold"]))

            # Fourth, convert to gromacs
            try:
                if 4 in steps:
                    if not is_automated:
                        input("Press enter to continue...")
                    chain_start, chain_end = convert_to_gromacs()
            except Exception as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [
                        index,
                        pdbid,
                        ligand_charge,
                        cofactor_charges,
                        ligand_count,
                        cofactor_count,
                        "Error",
                        "Error (Step 4)",
                    ]
                )
                continue

            print(colored("Step 4 completed successfully", "green", attrs=["bold"]))

            # Fifth, itp reindexing
            try:
                if 5 in steps:
                    if not is_automated:
                        input("Press enter to continue...")
                    itp_reindexing(chain_start, tuple(cofactors))
            except Exception as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [
                        index,
                        pdbid,
                        ligand_charge,
                        cofactor_charges,
                        ligand_count,
                        cofactor_count,
                        "Error",
                        "Error (Step 5)",
                    ]
                )
                continue

            print(colored("Step 5 completed successfully", "green", attrs=["bold"]))

            # Sixth, check in vmd
            try:
                if 6 in steps:
                    if not is_automated:
                        print(
                            colored(
                                f"Opening VMD for {index}.{pdbid}",
                                "yellow",
                                attrs=["bold"],
                            )
                        )
                        subprocess.run(
                            "vmd -e ../../scripts/state.vmd",
                            shell=True,
                            check=True,
                            executable=BASH,
                        )
            except Exception as error:
                print(colored(error, "red", attrs=["bold"]))
                writer.writerow(
                    [
                        index,
                        pdbid,
                        ligand_charge,
                        cofactor_charges,
                        ligand_count,
                        cofactor_count,
                        "Error",
                        "Error (Step 6)",
                    ]
                )
                continue

            print(colored("Step 6 completed successfully", "green", attrs=["bold"]))
            print(colored("All steps completed successfully", "green", attrs=["bold"]))
