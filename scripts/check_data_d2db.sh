#!/bin/bash 
d="/Shares/external-data/DAT_109__TCGA_unc_rnaseq_genefusion/by_tumor"
# branch1: cancer type
# branch2: rclbg vs nextgen2
# branch3: inside illumina (name)
# branch4: inside seqware-0.7.0_Mapsplice-0.7.4 (name_lane_barcode)

cd $d
for branch1 in ./*; do 
#for branch1 in HNSC; do
	if [ $branch1 != "._.DS_Store" ] && [ $branch1 != ".DS_Store" ]
	then
		path=$d/$branch1/datastore;
		echo "In "$path; cd $path
		for branch2 in nextgenout2 rclbg/nextgenout3; do
			if [ -d $branch2 ]; then
			
				path=$d/$branch1/datastore/$branch2/seqware-analysis/illumina;
				echo "In "$path; cd $path

				for branch3 in ./*; do

					path=$d/$branch1/datastore/$branch2/seqware-analysis/illumina/$branch3/seqware-0.7.0_Mapsplice-0.7.4;
					echo "In "$path; cd $path

					for branch4 in ./*; do

							echo $branch4 >> /Shares/work/DAT_109__TCGA_unc_rnaseq_genefusion/Data/d2db_files_present.log

					done # for_branch4

					#cd ../../
					cd $d/$branch1/datastore/$branch2/seqware-analysis/illumina

				done #for_branch3

				cd $d/$branch1/datastore
			fi # if branch2 exists

		done #for_branch2

		cd $d
		fi
done #for_branch_1

