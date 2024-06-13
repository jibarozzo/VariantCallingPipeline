#!/bin/bash
#SBATCH --job-name=rmdup\2-001
#SBATCH --mail-type=ALL # Valid values: BEGIN, END, FAIL, REQUEUE and ALL
#SBATCH --mail-user=cdong1@tulane.edu
#SBATCH --output=/lustre/project/kferris/Caroline_Dong/SHGL/2_demultiplex/rmdup.CMD2.L001.out
#SBATCH --error=/lustre/project/kferris/Caroline_Dong/SHGL/2_demultiplex/rmdup.CMD2.L001.error
#SBATCH --qos=normal
#SBATCH --time=12:00:00 		#: h:m:s max 24 hrs
#SBATCH --cpus-per-task=16   #: Cpus per Task
#SBATCH --nodes=1           #: Number of Nodes
#SBATCH --ntasks-per-node=1  #: Number of Tasks per Node
echo Start Job

### load modules
module load anaconda

### first step demultiplexing

python rmdup_molbarcodes.py -p CMD2-L001 -s fastq.gz

echo Finish CMD2 L001

module purge

echo End job