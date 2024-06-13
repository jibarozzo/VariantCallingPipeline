#!/bin/bash

#SBATCH --qos=normal
#SBATCH --job-name BAR34_ddRAD_rmdup
#SBATCH --error /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/errors/BAR34_ddRAD_rmdup.error     
#SBATCH --output /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/errors/BAR34_ddRAD_rmdup.output  
#SBATCH --time=23:00:00
#SBATCH --mem=128000 #Up to 256000, the maximum increases queue time
#SBATCH --nodes=1               # One is enough. When running MPIs, anywhere from 2-4 should be good.
#SBATCH --ntasks-per-node=1    # Number of Tasks per Node (MPI processes)
#SBATCH --cpus-per-task=20      # Number of cpu-cores (threads) per task (OMP threads)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=baponterolon@tulane.edu

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK # Set OMP_NUM_THREADS to the number of CPUs per task we asked for.

##Modules/Singularity
module load anaconda

# Basic session info
echo Start Job

echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

#Run Python script
python rmdup_molbarcodes.py -p "BAR34" -s fastq.gz

module purge

echo End Job