#!/bin/bash
#SBATCH --job-name=BAR_bam_diagnostics
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --output=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_1_diagnostics/output/bam_diagnostics_%A_%a.out
#SBATCH --error=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_1_diagnostics/errors/bam_diagnostics_%A_%a.err
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --mem=256000 #Up to 256000, the maximum increases queue time
#SBATCH --nodes=1            #: Number of Nodes
#SBATCH --ntasks-per-node=1  #: Number of Tasks per Node
#SBATCH --cpus-per-task=20   #: Number of threads per task
#SBATCH --array=1-474  # Job array (1-n) when n is number of unique samples that came off the sequencer. 499 total BUT 46 samples in BAR29 missing


### LOAD MODULES ###
module load bwa
module load samtools/1.10
module load java-openjdk

######################################################
### Diagnostics for Simple_pre_processing_workflow ###
######################################################

echo "Start Job"
cd /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/alignments_untrimmed

### ASSIGNING VARIABLES ###
P=$(find /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/alignments_untrimmed/* -type d \
    | sort \
    | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')

SAMPLE=$(echo $P | cut -d "/" -f 11) #Retrieves sample name
echo ${SAMPLE}


SEQID="bar_mim3" #project name and date for bam header
REF="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MimulusGuttatus_reference/MguttatusTOL_551_v5.0.fa"
THREADS=20
TMPDIR="/lustre/project/svanbael/TMPDIR" #designated storage folders for temporary files (should be empty at end)
PICARD="/share/apps/picard/2.20.7/picard.jar"


### BAM DIAGNOSTICS ###
### Using new Picard syntax
java -jar $PICARD ValidateSamFile \
      -I ${SAMPLE}/${SAMPLE}_markdup.bam \
      -VALIDATION_STRINGENCY STRICT \
      -REFERENCE_SEQUENCE $REF \
      -MODE SUMMARY \


module purge
echo "End Job"