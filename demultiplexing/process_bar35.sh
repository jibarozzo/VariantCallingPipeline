#!/bin/bash
#SBATCH --job-name=process35
#SBATCH --mail-type=ALL # Valid values: BEGIN, END, FAIL, REQUEUE and ALL
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --output=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/process_shortreads35.out
#SBATCH --error=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/process_shortreads35.error
#SBATCH --qos=normal
#SBATCH --time=24:00:00 		#: h:m:s max 24 hrs
#SBATCH --cpus-per-task=20   #: Cpus per Task
#SBATCH --nodes=1            #: Number of Nodes 18 max
#SBATCH --ntasks-per-node=1  #: Number of Tasks per Node

echo Start Job

### load modules
module load anaconda
module load stacks

### demultiplex by individual, join forward and reverse reads, uses Stacks

process_shortreads -1 filtered_forward.fastq -2 filtered_reverse.fastq -b barcodes_bar35.txt -o /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/BAR35/ -i fastq -y gzfastq -r -c -q -E phred33 --inline_null

module purge

echo End BAR35

echo End job