#....................USAGE......................#
# sh third_combineall.sh index pdbid

index=$1
pdbid=$2

#..................These paths to be modified by the users...........................#

work_path=${index}.${pdbid}
ff_location='~/anaconda3/envs/AmberTools21/dat/leap/cmd/'

cd $work_path/
cp tleap.in tleap_protoss.in

#............Changing the ligand naming in tleap file.........................#
sed -i "s/ligand_clean_h/ligand_clean_h_prot/g" ${work_path}/tleap_protoss.in
sed -i "s/ligand.frcmod/ligand_prot.frcmod/g"  ${work_path}/tleap_protoss.in

#............Changing the cofactor naming in tleap file.........................#
sed -i "s/cofactor_0\./cofactor_prot_0./g" ${work_path}/tleap_protoss.in
sed -i "s/cofactor_1\./cofactor_prot_1./g" ${work_path}/tleap_protoss.in
sed -i "s/cofactor_2\./cofactor_prot_2./g" ${work_path}/tleap_protoss.in

#............Changing the complex_solvated naming in tleap file.........................#
sed -i "s/complex_solvated\./complex_solvated_prot./g" ${work_path}/tleap_protoss.in

#............tleap header change.................................................#
first_line=$(sed -n 1p tleap_protoss.in)
replace_line="source ${ff_location}/leaprc.protein.ff14SB"
sed -i "s!${first_line}!${replace_line}!g" ${work_path}/tleap_protoss.in

second_line=$(sed -n 2p tleap_protoss.in)
replace_line="source ${ff_location}/leaprc.gaff2"
sed -i "s!${second_line}!${replace_line}!g" ${work_path}/tleap_protoss.in

third_line=$(sed -n 3p tleap_protoss.in)
replace_line="source ${ff_location}/leaprc.water.tip3p"
sed -i "s!${third_line}!${replace_line}!g" ${work_path}/tleap_protoss.in

tleap -s -f tleap_protoss.in
