#!/bin/bash
#SBATCH --job-name=BAR_bam_reassignment
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --output=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_2_reassignment/output/bam_reassigned_%A_%a.out
#SBATCH --error=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_2_reassignment/errors/bam_reassigned_%A_%a.err
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --mem=256000 #Up to 256000, the maximum increases queue time
#SBATCH --nodes=1            #: Number of Nodes
#SBATCH --ntasks-per-node=1  #: Number of Tasks per Node
#SBATCH --cpus-per-task=20   #: Number of threads per task
#SBATCH --array=1-475  # Job array (1-n) when n is number of unique samples that came off the sequencer. 499 total BUT 46 samples in BAR29 missing


### LOAD MODULES ###
module load java-openjdk

###################################################################
### Adding or replacing read group information in the bam files ###
###################################################################

echo "Start Job"
cd /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/alignments_untrimmed

### ASSIGNING VARIABLES ###
P=$(find /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/alignments_untrimmed/* -type d \
    | sort \
    | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')

SAMPLE=$(echo $P | cut -d "/" -f 11) #Retrieves sample name
HEADER=$(echo $P | cut -d "/" -f 11)
echo ${SAMPLE}


SEQID="bar_mim3" # Project name and date for bam header
REF="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MimulusGuttatus_reference/MguttatusTOL_551_v5.0.fa"
THREADS=20
TMPDIR="/lustre/project/svanbael/TMPDIR" # Designated storage folders for temporary files (should be empty at end)
PICARD="/share/apps/picard/2.20.7/picard.jar"

RGLB="lib1" # Library name (could be anything)
RGPL="ILLUMINA" # Sequencing platform
RGPU="UNKNOWN" # Generic identifier for the platform unit

# This tool enables the user to either replace all read groups in the input file or \
# to add a read group to each read that does not have a read group. Which could be the source of errors.
echo "Adding or replacing read group information"

 java -jar $PICARD AddOrReplaceReadGroups \
       I=${HEADER}/${SAMPLE}_markdup.bam \
       O=${HEADER}/${SAMPLE}_markdup_rrg.bam \
       RGID=${SEQID} \
       RGLB=$RGLB \
       RGPL=$RGPL \
       RGPU=$RGPU \
       RGSM= ${SAMPLE}
       VALIDATION_STRINGENCY=STRICT


module purge
echo "End Job"
