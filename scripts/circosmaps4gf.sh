# =======================================
# Author: Olga Nikolova
# E-mail: olga.nikolova@gmail.com
# =======================================
# Make gf circos maps
# =======================================
#!/bin/bash

# Read-in our MAF-like GF format
# Data and description/format available
# in Synapse: syn1896426
# Create input file for circos to plot
cd data4synapse
for file in *.tab; do 
#file in data4synapse/HNSC.tab; do 
  awk -F'\t' '{ print $8" " $10" "$10" "$9" "$11" "$11}' $file | sed 's/" c/c/g' | sed 's/"//g' > tmpfile; 
	awk -F'\t' '{ print $1}' $file | sort | uniq | wc -l > samplecount
	newfile=`basename $file .tab`
	cat "samplecount" "tmpfile" >$newfile.txt
done
cd ../

cp data4synapse/*.txt data4circos/.

# Run perl script which computes thickness
# and removes infrequent links

cd data4circos
for file in *.txt; do 
#for file in HNSC.circos.txt; do 
	`perl ../gf2bundles1.pl $file`; 
done
cd ../
cp data4circos/*circos.txt circos2/.

# Create circos maps
cd circos2
for file in *.circos.txt; do
#for file in 'HNSC.circos.txt'; do
	sed 's/file          = foo/file=\/Users\/olia\/Projects\/genefusion\/data4circos\/'$file'/' <circos.conf >tmp.conf
	/Users/olia/Local/circos-0.64/bin/circos -conf tmp.conf
	newfile=`basename $file .circos.txt`
	cp circos.png $newfile."circos.png"
	cp circos.svg $newfile."circos.svg"
done
