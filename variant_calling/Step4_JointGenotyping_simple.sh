#!/bin/bash
#SBATCH --job-name=GATK_pt2
#SBATCH --mail-type=ALL # Valid values: BEGIN, END, FAIL, REQUEUE and ALL.
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --output=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/output-joint-genotyping/joint_geno2_%A_%a.out
#SBATCH --error=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/error-joint-genotyping/joint_geno2_%A_%a.err
#SBATCH --qos=long
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=20   #: Cpus per Task
#SBATCH --nodes=1            #: Number of Nodes
#SBATCH --ntasks-per-node=1  #: Number of Tasks per Node

## create the error and output-joint-genotyping folders first 
## Do NOT create 2_GenomicsDB_data_store first. Error if this folder already exists

### LOAD MODULES ###
module load samtools/1.10
module load gatk/4.1.8.1
module load java-openjdk/1.8.0


#####################################################################
<<GATK_VariantCalling
This script uses the gatk toolkit version gatk/4.1.8.1.
Written in spring 2022 by Natalie Gonzalez

Modified by BAR 2024-06-13
Minor modifications:
- Updated samtools from version 1.10 to 1.16.1
- Assigned variables for the path to the reference genome and the temporary directory.

This is part 2 of the gatk variant calling process. This script
        calls the gatk GenomicDBImport tool and consolidates the VCF files
        produced in the previous script. It also takes a sample map file
        which is a tab-delimited text file containing the sample name
        followed by a tab and then followed by a path to the vcf file
        for that sample. An example would be:

        KGF_02  /path-to/4_NCBI_GATK_variant_calling/1_GVCF_files/KGF_02.g.vcf.gz

Once the samples have been  consolidated, GenotypeCVCFs is called. It can only be
                passed either 1 single sample vcf file, 1 multisample vcf file, or a
                GenomicsDB workspace created by GenomicsDBImport. This script uses the
                workspace option.

note: normally before begining this workflow we'd begin with a
        base quality recalibration step on the analysis ready bam
        files, and conclude with the gatk VariantRecalibrator, but
        both steps are being omitted because they require VCF files
        of already known variants for the populations being worked with
GATK_VariantCalling
#####################################################################

echo "Start Job"
echo "Start Consolidation"

### GLOBAL VARIABLES ###
WD="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/3_Genotyped_GVCFs" # Path to working directory
GDBWORKSPACE="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/2_GenomicsDB_data_store" # Path to GenomicsDB workspace
INTLIST="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/intervals.list" # Path to interval list
SAMPLEMAP="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/sample_name_map.txt" # Path to sample map file
TMPDIR="/lustre/project/svanbael/TMPDIR" # Path to temporary directory

BATCHSIZE=50 # Number of samples to consolidate at a time
THREADS=20 # Number of threads to use

OUT_NAME="BAR" # Output name for consolidated VCF

### CONSOLIDATING VCFS ###
#The GATK4 Best Practice Workflow for SNP and Indel calling uses GenomicsDBImport to merge GVCFs from multiple samples
#batchsize of 0 means no batches were used (i.e. readers for all samples will be opened at once)
#an interval list has to be passed to --L if whole genome is to be considered
gatk --java-options "-Xmx4g -Xms4g" \
	GenomicsDBImport \
	--genomicsdb-workspace-path $GDBWORKSPACE \
	--batch-size $BATCHSIZE \
	--L $INTLIST\
	--sample-name-map $SAMPLEMAP \
	--tmp-dir $TMPDIR \
	--reader-threads $THREADS

echo "End Consolidation"
#################################################################################################

echo "Start Joint Genotyping"

### JOINT GENOTYPING ###
# Move to the working directory

cd ${WD}

# GenotypeGVCFs is the tool to use when you have called variants separately
gatk  --java-options "-Xmx4g" GenotypeGVCFs\
    -R $REF\
    -V gendb:///lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/2_GenomicsDB_data_store/\
    -L $INTLIST\
    -O ${OUT_NAME}_joint_genotyping.vcf

echo "End Joint Genotyping"
#################################################################################################
echo "Printing Statistics"

### COUNTVARIANTS ###
# Counts variant records in VCF file regardless of filter status
gatk  --java-options "-Xmx4g" CountVariants\
    -V ${OUT_NAME}_joint_genotyping.vcf\
    --output ${OUT_NAME}_CountVariants

### VARIANTS TO TABLE ###
# Extracts fields from VCF to tab-delimited table
gatk  --java-options "-Xmx4g" VariantsToTable\
    -V ${OUT_NAME}_joint_genotyping.vcf\
    --output ${OUT_NAME}_VCF.table

echo "Finished Printing Statistics"

module purge
echo "End Job"
