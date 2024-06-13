#!/bin/bash
#SBATCH --job-name=HaplotypeCaller
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --output=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/output-variant-calling/var_call_%A_%a.out
#SBATCH --error=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/error-variant-calling/var_call_%A_%a.err
#SBATCH --qos=long
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1            #: Number of Nodes
#SBATCH --ntasks-per-node=1  #: Number of Tasks per Node
#SBATCH --cpus-per-task=20   #: Number of threads per task
#SBATCH --array=1-475   # Job array (1-n) when n is number of unique samples that came off the sequencer


### LOAD MODULES ###
module load samtools/1.10
module load gatk/4.1.8.1
module load java-openjdk/1.8.0

#####################################################################
<<GATK_VariantCalling_simple
This script uses the gatk toolkit version gatk/4.1.8.1 and
was written in spring 2022 by Natalie Gonzalez

This is part 1 of the gatk variant calling process. This script
        calls the gatk HaplotypeCaller in GVCF mode, producing a vcf
        file per sample containing genotype likelihoods for variant alleles.
        These files can be viewed in IGV. A TBI file is created for each VCF.

	This script is intended for samples that were not split across lanes and
	is therefore being referred to as the "simple" version. For samples split
	across lanes there's an alternative workflow.

note: normally before begining this workflow we'd begin with a
        base quality recalibration step on the analysis ready bam
        files, and conclude with the gatk VariantRecalibrator, but
        both steps are being omitted because they require VCF files
        of already known variants for the populations being worked with
GATK_VariantCalling_simple
#####################################################################

echo "Start Job"

### ASSIGNING VARIABLES ###
P=$(find /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/alignments_untrimmed/* -type d \
    | sort \
    | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')

# returns path: /lustre/project/kferris/Caroline_Dong/F2WY/3_preprocessing_TOLgenome/alignments_untrimmed/${SAMPLE}
SAMPLE=$(echo $P | cut -d "/" -f 11) #Retrieves sample name

### SETTING WORKING DIRECTORY WHERE VARIANT CALLING OUTPUTS WILL GO ###
cd /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/1_GVCF_files

#### HaplotypeCaller
gatk --java-options "-Xmx8g" HaplotypeCaller --tmp-dir /lustre/project/svanbael/TMPDIR \
   --reference /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MimulusGuttatus_reference/MguttatusTOL_551_v5.0.fa\
   --input /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/3_preprocessing/alignments_untrimmed/${SAMPLE}/${SAMPLE}_markdup.bam\
   --output ${SAMPLE}.g.vcf.gz\
   --emit-ref-confidence GVCF\
   --sample-name ${SAMPLE}


module purge

echo "End Job"
