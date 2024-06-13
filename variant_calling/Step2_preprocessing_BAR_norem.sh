#!/bin/bash
#SBATCH --job-name=BAR_preproces
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --output=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/preprocessingout/pre-processing_%A_%a.out
#SBATCH --error=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/preprocessingerr/pre-processing_%A_%a.err
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --mem=256000 #Up to 256000, the maximum increases queue time
#SBATCH --nodes=1            #: Number of Nodes
#SBATCH --ntasks-per-node=1  #: Number of Tasks per Node
#SBATCH --cpus-per-task=20   #: Number of threads per task
#SBATCH --array=1-2  # Job array (1-n) when n is number of unique samples that came off the sequencer. 499 total


### LOAD MODULES ###
module load bwa/0.7.17 
module load samtools/1.16.1
module load java-openjdk/1.8.0

#####################################################################
<<Simple_pre_processing_workflow
SIMPLE WORKFLOW- samples not split across lanes

Modified by BAR 2024-06-03
Major changes include:
- Optimizing the alignment workflow by piping the output of samtools `sort` directly into `fixmate`.
- Adding alternative alignment workflow with bwa, samtools and picard.
- Adding read group information using picard.jar::AddOrRepleaceReadGRoups as opossed to bwa mem in hopes of reducing errors.

Minor changes include:
- Using samtools v.1.16.1 as opposed to v.1.10
- bwa/0.7.17 as opposed to bwa
- java-openjdk/1.8.0 as opposed to java-openjdk


This script is designed to prepare samples for GATK varient calling. 
It begins with sequence files in seqdata.fq.gz or seqdata.fq format
The workflow is as follows:
perform the alignment with bwa -mem and assign readgroup information,
convert resulting .sam files to .bam files. The resulting file will
be an aligned paired end and sorted bam file.

Output file: ${SAMPLE}_aln_pe_fm_sorted.bam

NOTE: This script makes use of an arrray and by nature of a job array not all .bam
files will be generated at the same time. Make sure the entire job is done running
before moving forward. You can check on a job by using the "squeue" command in the 
terminal.
 
For second half- The workflow I followed uses the following documentation:
https://eriqande.github.io/eca-bioinf-handbook/alignment-of-sequence-data-to-a-reference-genome-and-associated-steps.html
In summary, the steps taken should be asigning read groups, aligning with bwa, and
processing with samtools


Input file:  ${SAMPLE}_aln_pe_sorted.bam
Output files:${SAMPLE}_merged.bam
                ${SAMPLE}_markdup.bam
                ${SAMPLE}_markdup.bam.bai
                ${SAMPLE}_markdup.bam.flagstat.txt

Simple_pre_processing_workflow
#####################################################################

echo "Start Job"
cd /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC

### ASSIGNING VARIABLES ###
# Array 1-n represents n number of samples and SLURM_ARRAY_TASK_ID represents elements in that array.
# The following line allows us to link the elements of the array with the FW and RV read files (R1/R2).
# The sort step should sort samples alphanumerically.
# Note: the find command searches recursively, grep command selects either read 1 or 2.
R1=$(find /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC \
    | grep 1.fq.gz \
    | sort \
    | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
R2=$(find /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC \
    | grep 2.fq.gz \
    | sort \
    | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')

# R1=$(find /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC/OPN_9_L_POOL.1.fq.gz \
#     | sort \
#     | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
# R2=$(find /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC/OPN_9_L_POOL.2.fq.gz \
#     | sort \
#     | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')


echo ${R1}
echo ${R2}
echo ${SLURM_ARRAY_TASK_ID}

# Sample prefix from the R1/R2 files, is in the sample_names.list (ex: KGF_02.._L3)
# This part changes based on the naming system. The user will have to modify this as needed
# The first "SAMPLE" command takes the file name in the xth field of the "/" delimiter, so you want to make sure that you have the full file path written and it is in that position. 
# So if you input "/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC/OPN_11_2_F2WY.2.fq.gz" you get "OPN_11_2_F2WY" for both sample and header in this code. 
SAMPLE=$(echo $R1 | cut -d "/" -f 10 | cut -d "." -f 1) # Retrieves first element before "."
HEADER=$(echo $R1 | cut -d "/" -f 10 | cut -d "." -f 1) # Retrieves first element before "." 

echo ${SAMPLE}
echo ${HEADER}

SEQID="bar_mim3" #project name and date for bam header
REF="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MimulusGuttatus_reference/MguttatusTOL_551_v5.0.fa"
THREADS=20
TMPDIR="/lustre/project/svanbael/TMPDIR" #designated storage folders for temporary files (should be empty at end)
PICARD="/share/apps/picard/2.20.7/picard.jar"

### Read group information
RGLB="lib1" # Library name (could be anything)
RGPL="ILLUMINA" # Sequencing platform
RGPU="UNKNOWN" # Generic identifier for the platform unit

#####################################################################
<<BWA_Alignment
This BWA aligmnent uses some of Caiti Heil's workflow from the
runSeqAlignVarientCall_20190130.sh script which I have stored
on my local machine under my desktop "Ferris Lab Materials" folder.
In an email conversation she mentioned that she did not trim the 
sequencing data because her lab got better alignments when they didn't trim.
I've combined some of her methods with the suggested workflow in:
https://eriqande.github.io/eca-bioinf-handbook/alignment-of-sequence-data-to-a-reference-genome-and-associated-steps.html
BWA_Alignment
#####################################################################

### SETTING WORKING DIRECTORY WHERE BWA OUTPUTS WILL GO ###
# Make alignments_untrimmed folder
cd /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/alignments_untrimmed
mkdir ${HEADER}  #makes a directory for each biological sample.

### RUNNING THE ALIGNMENT ON UNTRIMMED DATA ###
#fixmate -r flag removes unaligned reads
echo "Start Alignment"

#################################################################
### BWA alignment and SAMTOOLS use for read group information ###
#################################################################
#ORIGINAL CODE
#bwa mem -R '@RG\tID:'${SEQID}'\tSM:'${SAMPLE}'\tLB:lib1' -t ${THREADS} ${REF} ${R1} ${R2} \  # Aligning and adding read group information
#| samtools view -hb -@ ${THREADS} - \
#| samtools sort -n -T $TMPDIR -@ ${THREADS} - -o ${HEADER}/${SAMPLE}_aln_pe_sorted.bam \
#| samtools fixmate -rm -@ ${THREADS} ${HEADER}/${SAMPLE}_aln_pe_sorted.bam - \
#| samtools sort -T $TMPDIR -@ ${THREADS} - -o ${HEADER}/${SAMPLE}_aln_pe_fm_sorted.bam

# OPTIMIZED CODE -BAR 2024-05-22
# We avoid unnecessary writing and reading from disk by piping the output of samtools `sort` directly into `fixmate`.
# bwa mem -R '@RG\tID:'${SEQID}'\tSM:'${SAMPLE}'\tLB:lib1' -t ${THREADS} ${REF} ${R1} ${R2} \
# | samtools view -hb -@ ${THREADS} - \
# | samtools sort -n -T $TMPDIR -@ ${THREADS} - \
# | samtools fixmate -rm -@ ${THREADS} ${HEADER}/${SAMPLE}_aln_pe_sorted.bam - \
# | samtools sort -T $TMPDIR -@ ${THREADS} - -o ${HEADER}/${SAMPLE}_aln_pe_fm_sorted.bam


##################################################################################################
### Allignment with BWA, quality control with SAMTOOLS, and read group information with PICARD ###
##################################################################################################
# Adapted from @bergcollete's GitHub: YNP_GWAS/scripts/YNP4alignment.sh 


### Map reads to the genome
echo "Aligning ${SAMPLE} with bwa mem"
bwa mem -t $t ${REF} ${R1} ${R2} | \

### Quality filter and sort sam, making a bam file
echo "Making bam, quality filtering, and sorting for ${SAMPLE}"
samtools sort -n -T $TMPDIR -@ ${THREADS} - -o ${HEADER}/${SAMPLE}_aln_pe_sorted.bam | \
samtools fixmate -rm -@ ${THREADS} ${HEADER}/${SAMPLE}_aln_pe_sorted.bam - | \
samtools sort -T $TMPDIR -@ ${THREADS} - -o ${HEADER}/${SAMPLE}_aln_pe_fm_sorted.bam \

### Add read groups
echo "Adding read group information to ${SAMPLE}"
java -jar /home/thom_nelson/opt/picard.jar AddOrReplaceReadGroups \
    INPUT=${HEADER}/${SAMPLE}_aln_pe_fm_sorted.bam \
    OUTPUT=${HEADER}/${SAMPLE}_aln_pe_fm_rg_sorted.bam \
    RGSM=${rdgrp} \
    RGLB=NEBNext \
    RGPL=Illumina \
    RGPU=UNKNOWN \
    VALIDATION_STRINGENCY=LENIENT # adds read groups

echo "End Alignment"
##################################################################################################

### MARK AND REMOVE DUPLICATE READS ###
echo "Marking and removing duplicate reads"
java -jar $PICARD MarkDuplicates \
     INPUT=${HEADER}/${SAMPLE}_aln_pe_fm_sorted.bam \
     OUTPUT=${HEADER}/${SAMPLE}_markdup.bam \
     METRICS_FILE=${HEADER}/${SAMPLE}_dup_metrics.txt \
     ASSUME_SORTED=true \
     REMOVE_DUPLICATES=true \
     VALIDATION_STRINGENCY=LENIENT \
     TMP_DIR=$TMPDIR

echo "End Marking and removing duplicate reads"
### INDEXING THE BAM FILE & FLAG STATS###
echo "Indexing the bam file and calculating flag stats"
samtools index -@ ${THREADS} ${HEADER}/${SAMPLE}_markdup.bam

### LOOKING AT ALIGNMENT AGAIN ###
### `flagstat` counts the number of alignments for each FLAG type and calculates and prints statistics
samtools flagstat -@ ${THREADS} ${HEADER}/${SAMPLE}_markdup.bam \
   > ${HEADER}/${SAMPLE}_markdup.bam.flagstat.txt


module purge
echo "End Job"
