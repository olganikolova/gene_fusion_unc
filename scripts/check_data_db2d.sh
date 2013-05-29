#!/bin/bash 
while read name lane barcode title type; do
	[ "$name" == "name" ] && continue
	# create path to file
	# by_tumor/BLCA/130227_UNC9-SN296_0341_AC1RM4ACXX/seqware-0.7.0_Mapsplice-0.7.4/130227_UNC9-SN296_0341_AC1RM4ACXX_4_GTTTCG/fusion_junction.txt/
	bc=""
	if [ "$barcode" == "<NULL>" ]
	then
		bc=""
	else
		bc="_"$barcode

	fi

	file="by_tumor"/"$type"/"$name"_"$lane$bc"/"fusion_junction.txt"/"fusion_junction.txt"
	#echo $file

	if [ ! -f $file ]; then
		echo "File "$file" not found!" >> "db2d_files_missing.log"
	else
		echo $name"_"$lane$bc >> "db2d_files_present.log"
	fi
	#echo $type #"---"$name"---"$name"_"$lane$bc
	#echo $file
done < db.tab
