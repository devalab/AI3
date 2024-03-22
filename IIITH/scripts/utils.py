import os
import subprocess
from time import sleep

import numpy as np
import parmed as pmd
import requests
import selenium.webdriver
import selenium.webdriver.firefox.service
import wget
from prody import calcDistance, parsePDB, writePDB

BASH = "/bin/bash"
ZSH = "/bin/zsh"


def download_protoss_pdb(pdb_id, download_location):
    url = "https://proteins.plus/api/protoss_rest"
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
    }

    json_data = {
        "protoss": {
            "pdbCode": pdb_id,
        },
    }

    response = requests.post(
        url,
        headers=headers,
        json=json_data,
        timeout=60,
    )

    json_response = response.json()
    status_code = json_response["status_code"]

    if status_code not in (200, 202):
        raise Exception(f"Error: {json_response['message']}")

    location = json_response["location"]
    protoss_hash = location.split("/")[-1]
    get_url = url + "/" + protoss_hash
    response = requests.get(get_url, timeout=60)
    json_response = response.json()
    status_code = json_response["status_code"]

    while status_code != 200:
        sleep(5)
        response = requests.get(get_url, timeout=60)
        json_response = response.json()
        status_code = json_response["status_code"]
        if status_code not in (200, 202):
            raise Exception(f"Error: {json_response['message']}")

    pdb_url = json_response["protein"]
    print("--- Downloading PDB file from protoss ---")
    wget.download(pdb_url, out=f"{download_location}/{pdb_id}.pdb")
    print()
    print("--- PDB file downloaded ---")
    print()


def open_ligand_in_browser(ligand_id):
    url = f"https://www.rcsb.org/ligand/{ligand_id}"

    firefox_bin = "/snap/firefox/current/usr/lib/firefox/firefox"
    firefoxdriver_bin = "/snap/firefox/current/usr/lib/firefox/geckodriver"

    options = selenium.webdriver.firefox.options.Options()
    options.binary_location = firefox_bin

    service = selenium.webdriver.firefox.service.Service(
        executable_path=firefoxdriver_bin
    )

    browser = selenium.webdriver.Firefox(service=service, options=options)
    browser.get(url)


def parse_pdb(pdbfile):
    structure = parsePDB(pdbfile)
    name = structure.getResnames()[0]
    chain = structure.getChids()[0]
    res_idx = structure.getResnums()[0]

    return name, chain, res_idx


def diff(ls1, ls2):
    return list(list(set(ls1) - set(ls2)) + list(set(ls2) - set(ls1)))


def intersection(ls1, ls2):
    return list(set(ls1) & set(ls2))


def string_replace(filename, search, replace):
    # Read in the file
    with open(filename, "r") as file:
        filedata = file.read()

    # replace all occurrences of the required string
    filedata = filedata.replace(search, replace)

    # Write the file out again
    with open("file.txt", "w") as file:
        file.write(filedata)


def get_charge(pdbpath, script_path):
    subprocess.run(f"bash {script_path} {pdbpath}", shell=True, check=True, timeout=5)
    with open("charge.dat", "r") as f:
        for line in f:
            if line.startswith("Assigning partial charges"):
                info = line.split()
                for i in range(len(info)):
                    if info[i] == "charge":
                        charge = info[i + 1][:-1]
                        break
    print(f"# Info: Charge is: {charge}")
    subprocess.run("rm charge.dat", shell=True, check=True)
    return int(charge)


def get_atom_count(pdbfile):
    structure = parsePDB(pdbfile)
    return structure.numAtoms()


def pep_status():
    return os.path.exists("peptide.pdb")


def amino_status():
    return os.path.exists("aminoacid.pdb")


def cofact_status():
    return os.path.exists("cofactor_0.pdb")


def rename_str(mol, ch_id):
    # first_resid_name = {"ASP":"ASH","GLU":"GLH","LYS":"LYN","HIS":"HIE"}
    first_resid_name = {"ASP": "ASH", "GLU": "GLH", "HIS": "HIE"}

    if mol in ("peptide", "aminoacid"):
        name = f"{mol}_prot.pdb"
        pep_CYS_file = open("pep_CYS.dat", "w")
        structure = parsePDB(name)

    # Renaming the Histidine residues
    try:
        HIS = structure.select("resname HIS")
        HIS_resids = np.unique(HIS.getResnums())
        HD1 = list(np.unique(structure.select("resname HIS name HD1").getResnums()))
    except:
        print("No histidines with delta substitution")

    try:
        HE2 = list(np.unique(structure.select("resname HIS name HE2").getResnums()))
        HIP = intersection(HD1, HE2)
    except:
        print("No histidnes with epsilon substitution")

    try:
        for i in HIP:
            rename = structure.select("resnum {0:d}".format(i))
            rename.setResnames("HIP")
    except:
        print("No HIP")

    try:
        HD1 = list(np.unique(structure.select("resname HIS name HD1").getResnums()))
        for i in HD1:
            rename = structure.select("resnum {0:d}".format(i))
            rename.setResnames("HID")
    except:
        print("No HID")

    try:
        (structure.select("resname HIS")).setResnames("HIE")
    except:
        print("Residues HIE not exist")

    try:
        ASH = list(np.unique(structure.select("resname ASP name HD2").getResnums()))
        for i in ASH:
            rename = structure.select("resnum {0:d}".format(i))
            rename.setResnames("ASH")

    except:
        print("Residues ASH not exist")

    try:
        GLH = list(np.unique(structure.select("resname GLU name HE2").getResnums()))
        for i in GLH:
            rename = structure.select("resnum {0:d}".format(i))
            rename.setResnames("GLH")
    except:
        print("Residues GLH not exist")

    try:
        LYS = list(np.unique(structure.select("resname LYS").getResnums()))
        HZ1 = list(np.unique(structure.select("resname LYS name HZ1").getResnums()))
        LYN = diff(LYS, HZ1)
        print(LYN)
        for i in LYN:
            rename = structure.select("resnum {0:d}".format(i))
            rename.setResnames("LYN")
    except:
        print("Residues LYN not exist")

    try:
        CYS = list(np.unique(structure.select("resname CYS").getResnums()))
        for i in range(len(CYS)):
            for j in range(i + 1, len(CYS)):
                atom1 = structure.select("resnum {0:d} name SG".format(CYS[i]))
                atom2 = structure.select("resnum {0:d} name SG".format(CYS[j]))
                dist = calcDistance(atom1, atom2)
                if dist <= 2.90:
                    print("dist", dist)
                    rename = structure.select(
                        "resnum {0:d} {1:d}".format(CYS[i], CYS[j])
                    )
                    rename.setResnames("CYX")
                    print(CYS[i], CYS[j])

        SSBOND = list(structure.select("resname CYX name SG").getSerials())
        print("SSBOND", SSBOND)
        for i in range(len(SSBOND)):
            for j in range(i + 1, len(SSBOND)):
                atom1 = structure.select("serial {0:d}".format(SSBOND[i]))
                atom2 = structure.select("serial {0:d}".format(SSBOND[j]))
                dist = calcDistance(atom1, atom2)
                if dist <= 2.90:
                    R1 = structure.select("serial {0:d}".format(SSBOND[i])).getResnums()
                    R2 = structure.select("serial {0:d}".format(SSBOND[j])).getResnums()
                    print(R1, R2)
                    if name == "peptide_prot.pdb":
                        peptide_proto = (
                            "bond ligand_{0}.{1}.SG ligand_{0}.{2}.SG".format(
                                ch_id, int(R1), int(R2)
                            )
                        )
                    else:
                        peptide_proto = (
                            "bond receptor_{0}.{1}.SG receptor_{0}.{2}.SG".format(
                                ch_id, int(R1), int(R2)
                            )
                        )
                    print(peptide_proto)
                    pep_CYS_file.write("{0} \n".format(peptide_proto))

    except:
        print("Residues CYX not exist")

    # Handling the first residue of the peptide
    try:
        f_resname = structure.getResnames()[0]
        f_resid = structure.getResnums()[0]
        print(f_resname, f_resid)
        if f_resname in first_resid_name.keys():
            rename = structure.select("resnum {0:d}".format(f_resid))
            rename.setResnames(first_resid_name[f_resname])

        f_resname = structure.getResnames()[-1]
        f_resid = structure.getResnums()[-1]
        print(f_resname, f_resid)
        if f_resname in first_resid_name.keys():
            print(first_resid_name.keys())
            rename = structure.select("resnum {0:d}".format(f_resid))
            rename.setResnames(first_resid_name[f_resname])
    except:
        print("Either first residue or last residue do not require renaming")

    writePDB(name, structure.select("noh"))


def file_map(old_file, newfile, mol, file1):
    d_file = parsePDB(file1)
    if mol == "lig" or mol == "aminoacid":
        s = parsePDB(old_file)
        ch_id = s.getChids()[0]
        res_id = s.getResnums()[0]
        res_name = s.getResnames()[0]
        if mol == "aminoacid":
            sel = d_file.select(
                "chain {0} and resid {1} and resname {2}".format(
                    ch_id, res_id, res_name
                )
            )
        elif mol == "lig":
            sel = d_file.select("chain {0} and resid {1}".format(ch_id, res_id))

        writePDB(newfile, sel)
        os.system('sed -i "s/\'/Z/g" {}'.format(newfile))

        print("Name of the ligand", res_name)
        if mol == "aminoacid":
            rename_str(mol, ch_id)
    elif mol == "peptide":
        s = parsePDB(old_file)
        ch_id = s.getChids()[0]
        sel = d_file.select("chain {0} and protein".format(ch_id))
        writePDB(newfile, sel)
        rename_str(mol, ch_id)


def minimization_check():
    output = subprocess.check_output(
        "tail -n 20 em.log", shell=True, encoding="utf-8"
    ).split("\n")
    potential_energy = 0
    for out in output:
        if "Potential Energy" in out:
            potential_energy = float(out.split("=")[1].strip())
    if potential_energy >= 0:
        raise IndexError
    return potential_energy


def amber_to_gmx(cofactors):
    os.system("rm -rf complex_solvated_prot.top complex_solvated_prot.gro")
    amber = pmd.load_file("complex_solvated_prot.prmtop", "complex_solvated_prot.pdb")
    amber.save("complex_solvated_prot.top", format="gromacs")
    amber.save("complex_solvated_prot.gro")

    structure = parsePDB("complex_solvated_prot.pdb")
    chain_end_sel = structure.select("protein and name OXT")
    chain_end = chain_end_sel.getIndices() + 1
    num_chains = len(chain_end)
    chain_start = [1] * num_chains
    for i in range(1, num_chains):
        chain_start[i] = chain_end[i - 1] + 1

    if cofact_status() and amino_status():
        cof_select_last = structure.select(f"resname {cofactors[-1]}")
        chain_start[-1] = cof_select_last.getIndices()[-1] + 2

    print(chain_start)
    print(chain_end)

    os.system("echo q | gmx make_ndx -f complex_solvated_prot.gro -o indx.ndx")

    for i in range(num_chains):
        if not pep_status() and not amino_status():
            chain_id = i + 1
            print(
                "Now will run: gmx make_ndx -f complex_solvated_prot.gro -n indx.ndx -o indx.ndx"
            )
            os.system(
                'printf "%s\n" "a {0}-{1}  & ! a H* & 5" "q" | gmx make_ndx -f complex_solvated_prot.gro -n indx.ndx -o indx.ndx'.format(
                    chain_start[i], chain_end[i]
                )
            )
            group_id = int(
                subprocess.check_output('grep "\[" indx.ndx | wc -l', shell=True)
            )

            os.system(
                'printf "%s\n" "name {0} chain_{1}" "q" |gmx make_ndx -f complex_solvated_prot.gro -n indx.ndx -o indx.ndx'.format(
                    group_id - 1, chain_id
                )
            )
            print(
                "Now will run: gmx genrestr -f complex_solvated_prot.gro -n indx.ndx -o posre_chain_{}.itp -fc 4000 4000 4000".format(
                    chain_id
                )
            )
            os.system(
                'echo "{0}" |gmx genrestr -f complex_solvated_prot.gro -n indx.ndx -o posre_chain_{1}.itp -fc 4000 4000 4000'.format(
                    group_id - 1, chain_id
                )
            )
        else:
            if i != num_chains - 1:
                chain_id = i + 1
                print(
                    "Now will run: gmx make_ndx -f complex_solvated_prot.gro -n indx.ndx -o indx.ndx"
                )
                os.system(
                    'printf "%s\n" "a {0}-{1}  & ! a H* & 5" "q" | gmx make_ndx -f complex_solvated_prot.gro -n indx.ndx -o indx.ndx'.format(
                        chain_start[i], chain_end[i]
                    )
                )
                group_id = int(
                    subprocess.check_output('grep "\[" indx.ndx | wc -l', shell=True)
                )

                os.system(
                    'printf "%s\n" "name {0} chain_{1}" "q" |gmx make_ndx -f complex_solvated_prot.gro -n indx.ndx -o indx.ndx'.format(
                        group_id - 1, chain_id
                    )
                )
                print(
                    "Now will run: gmx genrestr -f complex_solvated_prot.gro -n indx.ndx -o posre_chain_{}.itp -fc 4000 4000 4000".format(
                        chain_id
                    )
                )
                os.system(
                    'echo "{0}" |gmx genrestr -f complex_solvated_prot.gro -n indx.ndx -o posre_chain_{1}.itp -fc 4000 4000 4000'.format(
                        group_id - 1, chain_id
                    )
                )
            else:
                chain_id = i + 1
                os.system(
                    'printf "%s\n" "a {0}-{1}  & ! a H* & 5" "q" | gmx make_ndx -f complex_solvated_prot.gro -n indx.ndx -o indx.ndx'.format(
                        chain_start[i], chain_end[i]
                    )
                )
                group_id = int(
                    subprocess.check_output('grep "\[" indx.ndx | wc -l', shell=True)
                )
                os.system(
                    'printf "%s\n" "name {0} chain_{1}" "q" |gmx make_ndx -f complex_solvated_prot.gro -n indx.ndx -o indx.ndx'.format(
                        group_id - 1, chain_id
                    )
                )
                os.system(
                    'echo "{0}" |gmx genrestr -f complex_solvated_prot.gro -n indx.ndx -o posre_ligand_prot.itp -fc 4000 4000 4000'.format(
                        group_id - 1
                    )
                )
                os.system("cp posre_ligand_prot.itp temp_lig.itp")

    return tuple(chain_start), tuple(chain_end)


def itp_reindexing(chain_start, cof):
    f1 = open("complex_solvated_prot.top", "r")
    fl1 = f1.readlines()
    f1.close()
    num_chains = len(chain_start)
    num_cofactors = len(cof)
    for chain_id in range(2, num_chains + 1):
        print("working on the itp file of chain {}.\n".format(chain_id))
        if not pep_status() and not amino_status():
            print("posre_chain_{}.itp".format(chain_id))
            f2 = open("posre_chain_{}.itp".format(chain_id), "r")
        else:
            if chain_id != num_chains:
                print("posre_chain_{}.itp".format(chain_id))
                f2 = open("posre_chain_{}.itp".format(chain_id), "r")
            else:
                print("Working on ligand posre")
                f2 = open("posre_ligand_prot.itp", "r")

        fl2 = f2.readlines()
        f = open("temp.itp", "w")
        n = chain_start[chain_id - 1] - 1

        for i in range(len(fl2)):
            try:
                temp = int(fl2[i].split()[0]) - int(n)
                f.write(
                    "{0:6d} {1:5d} {2:11d} {3:11d} {4:11d}\n".format(
                        temp,
                        int(fl2[i].split()[1]),
                        int(fl2[i].split()[2]),
                        int(fl2[i].split()[3]),
                        int(fl2[i].split()[4]),
                    )
                )

            except:
                f.write(fl2[i])

        f2.close()
        f.close()

        if not pep_status() and not amino_status():
            os.system("mv temp.itp posre_chain_{}.itp".format(chain_id))
            os.system("cp ligand_prot.acpype/posre_ligand_prot.itp .")
        else:
            if chain_id != num_chains:
                os.system("mv temp.itp posre_chain_{}.itp".format(chain_id))
            else:
                os.system("mv temp.itp posre_ligand_prot.itp")

    system = []

    for i, line in enumerate(fl1):
        if "moleculetype" in line:
            system.append(fl1[i + 2].split()[0])

    print(system)
    system_count = []

    j = 0
    for i in range(len(system)):
        system_count.append(int(fl1[-1 + j].split()[1]))
        j = j - 1

    f2 = open("temp.top", "w")

    chain_count = 0
    cof_flag = []

    for i in range(num_cofactors):
        cof_flag.append(0)

    lig_flag = 0

    if not pep_status() and not amino_status():
        protein_itp = num_chains
        LIG_string = "LIG"
    elif pep_status():
        protein_itp = num_chains - 1
        LIG_string = "system{}".format(num_chains)
    elif amino_status():
        protein_itp = num_chains - 1
        st = parsePDB("aminoacid_prot.pdb")
        LIG_string = st.getResnames()[0]

    print("Number of chains:", num_chains, "LIG_string:", LIG_string)

    for i in range(len(system)):
        if system[i] == LIG_string:
            next_to_lig = system[i + 1]
            print(next_to_lig)

    protein_posre_string = ""
    skip_itp = [0] * len(system)

    for i in range(len(fl1)):
        if (
            "moleculetype" in fl1[i]
            and "system" in fl1[i + 2]
            and chain_count > 0
            and chain_count <= protein_itp
            and "#ifdef" not in fl1[i - 4]
        ):
            f2.write("#ifdef POSRES_{}\n".format(chain_count))
            f2.write('#include "posre_chain_{}.itp"\n'.format(chain_count))
            f2.write("#endif\n\n")
            protein_posre_string = (
                protein_posre_string + "-DPOSRES_{}".format(chain_count) + " "
            )
            chain_count = chain_count + 1

        if int(system_count[chain_count - 1]) > 1:
            if not pep_status() and not amino_status():
                if int(skip_itp[chain_count - 1]) == 0:
                    if (
                        "moleculetype" in fl1[i]
                        and system[chain_count] in fl1[i + 2]
                        and chain_count > 0
                        and chain_count <= protein_itp
                        and "#ifdef" not in fl1[i - 4]
                    ):
                        f2.write("#ifdef POSRES_{}\n".format(chain_count))
                        f2.write('#include "posre_chain_{}.itp"\n'.format(chain_count))
                        f2.write("#endif\n\n")
                        protein_posre_string = (
                            protein_posre_string
                            + "-DPOSRES_{}".format(chain_count)
                            + " "
                        )
                        chain_count = chain_count + 1
                        skip_itp[chain_count - 1] = skip_itp[chain_count - 1] + 1

        if "moleculetype" in fl1[i] and "system" in fl1[i + 2] and chain_count == 0:
            chain_count = chain_count + 1

        if num_cofactors > 0:
            if (
                "moleculetype" in fl1[i]
                and cof[0] in fl1[i + 2]
                and chain_count == protein_itp
            ):
                f2.write("#ifdef POSRES_{}\n".format(chain_count))
                f2.write('#include "posre_chain_{}.itp"\n'.format(chain_count))
                f2.write("#endif\n\n")
                protein_posre_string = (
                    protein_posre_string + "-DPOSRES_{}".format(chain_count) + " "
                )
                chain_count = chain_count + 1
        else:
            if (
                "moleculetype" in fl1[i]
                and LIG_string in fl1[i + 2]
                and chain_count == protein_itp
            ):
                f2.write("#ifdef POSRES_{}\n".format(chain_count))
                f2.write('#include "posre_chain_{}.itp"\n'.format(chain_count))
                f2.write("#endif\n\n")
                protein_posre_string = (
                    protein_posre_string + "-DPOSRES_{}".format(chain_count) + " "
                )
                chain_count = chain_count + 1

        if num_cofactors > 0:
            for cof_ind in range(num_cofactors):
                if "moleculetype" in fl1[i] and cof[cof_ind] in fl1[i + 2]:
                    cof_flag[cof_ind] = 1
                    print("Working on cofactor {}".format(cof_ind))
                try:
                    if (
                        "moleculetype" in fl1[i]
                        and cof[cof_ind + 1] in fl1[i + 2]
                        and cof_flag[cof_ind] == 1
                    ):
                        f2.write("#ifdef POSRES_COF{}\n".format(cof_ind))
                        f2.write(
                            '#include "posre_cofactor_prot_{}.itp"\n'.format(cof_ind)
                        )
                        f2.write("#endif\n\n")
                except IndexError:
                    if (
                        "moleculetype" in fl1[i]
                        and LIG_string in fl1[i + 2]
                        and cof_flag[cof_ind] == 1
                    ):
                        f2.write("#ifdef POSRES_COF{}\n".format(num_cofactors - 1))
                        f2.write(
                            '#include "posre_cofactor_prot_{}.itp"\n'.format(
                                num_cofactors - 1
                            )
                        )
                        f2.write("#endif\n\n")

        if "moleculetype" in fl1[i] and LIG_string in fl1[i + 2]:
            lig_flag = 1

        if (
            "moleculetype" in fl1[i]
            and lig_flag == 1
            and fl1[i + 7].split()[3] == next_to_lig
        ):
            print(fl1[i + 7])
            f2.write("#ifdef POSRES_LIG\n")
            f2.write('#include "posre_ligand_prot.itp"\n')
            f2.write("#endif\n\n")
            f2.write(fl1[i])
        else:
            f2.write(fl1[i])

    f2.close()

    replace_string = protein_posre_string
    search_string = "protein_posres"
    string_replace("em.mdp", search_string, replace_string)
    string_replace("nvt.mdp", search_string, replace_string)

    os.system("mv temp.top complex_solvated_prot.top")

    os.system(
        "gmx grompp -f em.mdp -o em.tpr -c complex_solvated_prot.gro -r complex_solvated_prot.gro -p complex_solvated_prot.top -n indx.ndx -maxwarn 1"
    )
    os.system("gmx mdrun -deffnm em -v")
    minimization_check()
