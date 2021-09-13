#!/bin/bash

folder_out=1
END=1000000 # ps
#END=60000 # ps
step=5000 # ps
start1=50000 #ps

dir="md5ns_pdbs_CA_out"

mkdir ${dir}

for ((i=$start1;i<$END;i=i+$step)); # do not use END frame 
do

	echo $i
	echo 19 3 | gmx trjconv -s ./md5ns_frame50000.tpr -f ../md1us1_dt100.trr -o ${dir}/frame_${i}ps_CA_no_pbc.xtc -center -pbc nojump -ur compact -b ${i} -e ${i} -n after_dex_removal.ndx #>&/dev/null

	echo 10 3 | gmx trjconv -s ./out_xtc_md1/ca_frames/md5ns_frame50000_CA.tpr -f ${dir}/frame_${i}ps_CA_no_pbc.xtc -o ${dir}/frame_${i}ps_CA_fit.pdb -fit rot+trans -n ./out_xtc_md1/ca_frames/index.ndx #>&/dev/null


done

