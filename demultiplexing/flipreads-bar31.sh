#!/bin/bash
#SBATCH --job-name=flip31
#SBATCH --mail-type=ALL # Valid values: BEGIN, END, FAIL, REQUEUE and ALL
#SBATCH --mail-user=cdong1@tulane.edu
#SBATCH --output=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/flipBAR31.out
#SBATCH --error=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/flipBAR31.error
#SBATCH --qos=normal
#SBATCH --time=24:00:00 		#: h:m:s max 24 hrs
#SBATCH --cpus-per-task=20   #: Cpus per Task
#SBATCH --nodes=1            #: Number of Nodes 18 max
#SBATCH --ntasks-per-node=1  #: Number of Tasks per Node

echo Start Job

### load modules
module load anaconda

### first step demultiplexing

python2 flip2BeRAD.py -c TGCAG -f BAR31.rmdup.1.fastq -r BAR31.rmdup.2.fastq -b ddrad_barcode_seqs.txt -m 0 -o 2

## Move output files into another folder for further demultiplexing by individual
mv -i barcode_no_cut_reverse.fastq BAR31/
mv -i barcode_no_cut_forward.fastq BAR31/
mv -i filtered_forward.fastq BAR31/
mv -i filtered_reverse.fastq BAR31/
mv -i nobarcodes_reverse.fastq BAR31/
mv -i nobarcodes_forward.fastq BAR31/
mv -i log.txt BAR31/

module purge

echo End BAR31

echo End job