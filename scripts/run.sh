#!/bin/bash

# =======================================
# Author: Olga Nikolova
# E-mail: olga.nikolova@gmail.com
# =======================================
# Make gene fusion data from unc actionable
# =======================================

# Step 1: Move all files into one dir per cancer
d=`pwd`
mcdir $d/mapped2TCGA_annot
cd by_tumor
for ct in *; do
    for subdir in nextgenout2 rclbg/nextgenout3; do
        cd $d/by_tumor/$ct/datastore/$subdir/seqware-analysis/illumina
        for name in *; do
        if [ ! -d $d/mapped2TCGA_annot/$ct ]
        then
            mkdir $d/mapped2TCGA_annot/$ct
        fi
        mv $d/by_tumor/$ct/datastore/$subdir/seqware-analysis/illumina/$name/seqware-0.7.0_Mapsplice-0.7.4/* $d/mapped2TCGA_annot/$ct/.
        rm -r $d/by_tumor/$ct/datastore
        done
    done
done

# Step 2: Format files
while read name lane barcode title type; do
    [ "$name" == "name" ] && continue

    # create path to file
    bc=""
    if [ "$barcode" == "<NULL>" ]
    then
        bc=""
    else
        bc="_"$barcode
    fi

    file="by_tumor"/"$type"/"$name"_"$lane$bc"/"fusion_junction.txt"/"fusion_junction.txt"

    # Copy the mapped existing files into mapped2TCGA by the unique filename
    if [ -f $file ]; then
        if [ ! -d mapped2TCGA/$type ]
        then
            mkdir mapped2TCGA/$type
        fi

        # Create new filename: $name_$lane$bc.txt
        newfile=mapped2TCGA/$type/$name"_"$lane$bc".txt"
        cp $file $newfile

        if [ ! -d mapped2TCGA_annot/$type ]
        then
            mkdir mapped2TCGA_annot/$type
        fi

        # 3. Add the each files columns:
        # title: TCGA barcode
        # type
        # name
        # lane
        # barcode
        # and store in mapped2TCGA_annot/$type
        awk -v newcolumns=$title"\t"$type"\t"$name"\t"$lane"\t"$barcode"\t" '{print newcolumns,$0}' <$newfile >mapped2TCGA_annot/$type/$name"_"$lane$bc"_annot.txt"
    fi

done < db.tab

# Step 3: Cat all files within tumor type
cd mapped2TCGA_annot
for type in ./*;do
    cat $type/*_annot.txt > ../mapped2TCGA_final/$type.txt
done
cd ../

# Step 3**: Split the ChrA~ChrB column into donorChromosome and acceptorChromosome
cd mapped2TCGA_final
for f in ./*; do
    sed -i 's/~/\t/g' $f
done
cd ../

