#!/bin/bash

Red='\e[1;31m'
Yel='\e[1;33m'
RCol='\e[0m'


read -p "$(echo ${Yel} index: `echo '\n> '` $RCol)" index
read -p "$(echo ${Yel} pdbid: `echo '\n> '` $RCol)" pdbid
read -p "$(echo ${Yel} N_cof: `echo '\n> '` $RCol)" N_cof


# ligand parameter genreation
lig_file_name='ligand_prot.pdb'
if [ -e ${lig_file_name} ]; then
	read -p "$(echo ${Yel} ligand_charge: `echo '\n> '` $RCol)" ligand_charge
	antechamber -fi pdb -fo prepi -i ligand_prot.pdb -o ligand_clean_h_prot.prepi -rn LIG -nc $ligand_charge -c bcc -pf y -at gaff2  -ek "maxcyc=0"
	if [ ! -e ligand_clean_h_prot.prepi ]
	then
		antechamber -fi pdb -fo prepi -i ligand_prot.pdb -o ligand_clean_h_prot.prepi -rn LIG -nc $ligand_charge -c bcc -pf y -at gaff2  
		antechamber -fi pdb -fo mol2 -i ligand_prot.pdb -o ligand_clean_h_prot.mol2 -rn LIG -c bcc -nc $ligand_charge -pf y -at gaff2   
		parmchk2 -f prepi -i ligand_clean_h_prot.prepi -o ligand_prot.frcmod
	elif [ ! -e ligand_clean_h_prot.prepi ]
	then 
                antechamber -fi pdb -fo prepi -i ligand_prot.pdb -o ligand_clean_h_prot.prepi -rn LIG -nc $ligand_charge -c bcc -pf y -at gaff2 -ek "maxcyc=0 ndiis_attempts=700"
                antechamber -fi pdb -fo mol2 -i ligand_prot.pdb -o ligand_clean_h_prot.mol2 -rn LIG -c bcc -nc $ligand_charge -pf y -at gaff2 -ek "maxcyc=0 ndiis_attempts=700"
                parmchk2 -f prepi -i ligand_clean_h_prot.prepi -o ligand_prot.frcmod
	else
		antechamber -fi pdb -fo mol2 -i ligand_prot.pdb -o ligand_clean_h_prot.mol2 -rn LIG -c bcc -nc $ligand_charge -pf y -at gaff2  -ek "maxcyc=0"
		parmchk2 -f prepi -i ligand_clean_h_prot.prepi -o ligand_prot.frcmod
	fi

	acpype -i ligand_clean_h_prot.mol2 -c user -n $ligand_charge -a gaff2 -b ligand_prot -o gmx
	cp ligand_prot.acpype/posre_ligand_prot.itp .
	#rm sqm.in sqm.out sqm.pdb
else
	echo "${Red} ligand is not an organic molecule ${RCol}"
fi

# cofactor parameter generation
N_cof=$(ls -l cofactor_prot_*.pdb| wc -l )
if [ ${N_cof} > 0 ]; then
	i=0
	while [ "$i" -lt "$N_cof" ]
	do
		echo "Working on cofactor ${i}"
		read -p "$(echo ${Yel} cofactor_charge: `echo '\n:'` $RCol)" cofactor_charge
		antechamber -fi pdb -fo prepi -i cofactor_prot_${i}.pdb -o cofactor_prot_${i}.prepi -c bcc -pf y -nc $cofactor_charge -at gaff2 -ek "maxcyc=0"
		if [ ! -e cofactor_prot_${i}.prepi ]
		then
			antechamber -fi pdb -fo prepi -i cofactor_prot_${i}.pdb -o cofactor_prot_${i}.prepi -c bcc -pf y -nc $cofactor_charge -at gaff2 
			antechamber -fi pdb -fo mol2 -i cofactor_prot_${i}.pdb -o cofactor_prot_${i}.mol2 -c bcc -pf y -nc $cofactor_charge -at gaff2 
			parmchk2 -f prepi -i cofactor_prot_${i}.prepi -o cofactor_prot_${i}.frcmod
			rm sqm.in sqm.out sqm.pdb
		elif [ ! -e cofactor_prot_${i}.prepi ]
		then
                        antechamber -fi pdb -fo prepi -i cofactor_prot_${i}.pdb -o cofactor_prot_${i}.prepi -c bcc -pf y -nc $cofactor_charge -at gaff2 -ek "maxcyc=0 ndiis_attempts=700"
                        antechamber -fi pdb -fo mol2 -i cofactor_prot_${i}.pdb -o cofactor_prot_${i}.mol2 -c bcc -pf y -nc $cofactor_charge -at gaff2 -ek "maxcyc=0 ndiis_attempts=700"
                        parmchk2 -f prepi -i cofactor_prot_${i}.prepi -o cofactor_prot_${i}.frcmod

		else
			antechamber -fi pdb -fo mol2 -i cofactor_prot_${i}.pdb -o cofactor_prot_${i}.mol2 -c bcc -pf y -nc $cofactor_charge -at gaff2 -ek "maxcyc=0"
			parmchk2 -f prepi -i cofactor_prot_${i}.prepi -o cofactor_prot_${i}.frcmod
			rm sqm.in sqm.out sqm.pdb
		fi

		acpype -i cofactor_prot_${i}.mol2 -c user -n $cofactor_charge -a gaff2 -b cofactor_prot_${i} -o gmx
		cp cofactor_prot_${i}.acpype/posre_cofactor_prot_${i}.itp .


		i=$(($i+1))
	done

else
	echo "${Red} Cofactor doesn't exist ${RCol}"
fi







