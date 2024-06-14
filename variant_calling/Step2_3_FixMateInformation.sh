#!/bin/bash
#SBATCH --job-name=BAR_bam_fixedmate
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --output=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_2_bam_fixedmate/output/bam_fixedmate_%A_%a.out
#SBATCH --error=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_2_fixedmate/errors/bam_fixedmate_%A_%a.err
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --mem=256000 #Up to 256000, the maximum increases queue time
#SBATCH --nodes=1            #: Number of Nodes
#SBATCH --ntasks-per-node=1  #: Number of Tasks per Node
#SBATCH --cpus-per-task=20   #: Number of threads per task
#SBATCH --array=1-475  # Job array (1-n) when n is number of unique samples that came off the sequencer. 499 total BUT 46 samples in BAR29 missing


### LOAD MODULES ###
module load java-openjdk

###############################################
### Fixin mate information in the bam files ###
###############################################

echo "Start Job"

### NAVIGATING TO THE DIRECTORY CONTAINING THE BAM FILES ###
### WD should be the directory containing the bam files ###
WD="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/alignments_untrimmed/" 
cd ${WD}

### ASSIGNING VARIABLES ###
P=$(find ${WD}* -type d \
    | sort \
    | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')

SAMPLE=$(echo $P | cut -d "/" -f 11) #Retrieves sample name

echo ${SAMPLE}

SEQID="bar_mim3" # Project name and date for bam header
REF="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MimulusGuttatus_reference/MguttatusTOL_551_v5.0.fa" # Path to reference genome
THREADS=20 # Number of threads to use
TMPDIR="/lustre/project/svanbael/TMPDIR" # Designated storage folders for temporary files (should be empty at end)
PICARD="/share/apps/picard/2.20.7/picard.jar" # Path to picard

### FIXING MATE INFORMATION IN BAM FILES ###
echo "Fixing mate information in ${SAMPLE} bam file"

java -jar $PICARD FixMateInformation \
       -I ${SAMPLE}/${SAMPLE}_markdup_rrg.bam \
       -O ${SAMPLE}/${SAMPLE}_markdup_rrg_fm.bam \
       -ADD_MATE_CIGAR true

echo "Finished fixing mate information in ${SAMPLE} bam file"

module purge
echo "End Job"