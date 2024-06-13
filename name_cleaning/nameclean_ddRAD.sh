#!/bin/bash
# for file in *.fastq.gz
# do
#  echo "Unziping"
#  gzip -d  "$file"
# done

#Rename the files

# Iterate over all files with the pattern *R1_001.fastq.gz
for file in *R1_001.fastq.gz; do
            echo "Renaming R1s"
            # Extract the prefix before the first underscore and append .1.fastq.gz
            newname=$(echo "$file" | cut -d_ -f1).1.fastq.gz
            echo "Renaming $file as $newname"
            # Uncomment the line below to actually perform the renaming
            mv "$file" "$newname"
done


# Iterate over all files with the pattern *R3_001.fastq.gz
for file in *R3_001.fastq.gz; do
            echo "Renaming R2s"
            # Extract the prefix before the first underscore and append .2.fastq.gz
            newname=$(echo "$file" | cut -d_ -f1).2.fastq.gz
            echo "Renaming $file as $newname"
            # Uncomment the line below to actually perform the renaming
            mv "$file" "$newname"
done
    
    
# Iterate over all files with the pattern *I1_001.fastq.gz
for file in *I1_001.fastq.gz; do
            echo "Renaming I7s"
            # Extract the prefix before the first underscore and append .I7.fastq.gz
            newname=$(echo "$file" | cut -d_ -f1).i7.fastq.gz
            echo "Renaming $file as $newname"
            # Uncomment the line below to actually perform the renaming
            mv "$file" "$newname"
done


# Iterate over all files with the pattern *R2_001.fastq.gz
for file in *R2_001.fastq.gz; do
            echo "Renaming I5s"
            # Extract the prefix before the first underscore and append .I7.fastq.gz
            newname=$(echo "$file" | cut -d_ -f1).i5.fastq.gz
            echo "Renaming $file as $newname"
            # Uncomment the line below to actually perform the renaming
            mv "$file" "$newname"
done

# Code augmented with ChatGPT's help.
