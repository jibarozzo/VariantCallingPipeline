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
This script uses the gatk toolkit version gatk/4.1.8.1 and
was written in spring 2021 by Natalie Gonzalez

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

### CONSOLIDATING VCFS ###
#The GATK4 Best Practice Workflow for SNP and Indel calling uses GenomicsDBImport to merge GVCFs from multiple samples
#batchsize of 0 means no batches were used (i.e. readers for all samples will be opened at once)
#an interval list has to be passed to --L if whole genome is to be considered
gatk --java-options "-Xmx4g -Xms4g" \
	GenomicsDBImport \
	--genomicsdb-workspace-path /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/2_GenomicsDB_data_store \
	--batch-size 50 \
	--L /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/intervals.list \
	--sample-name-map /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/sample_name_map.txt \
	--tmp-dir /lustre/project/svanbael/TMPDIR \
	--reader-threads 20

#echo "End Consolidation"
#################################################################################################

echo "Start Joint Genotyping"

### JOINT GENOTYPING ###
#change working directory
#mkdir /lustre/project/kferris/Diana_Tataru/F4_F5/4_F4_F5GATK_variant_calling/3_Genotyped_GVCFs
cd /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/3_Genotyped_GVCFs

gatk  --java-options "-Xmx4g" GenotypeGVCFs\
    -R /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MimulusGuttatus_reference/MguttatusTOL_551_v5.0.fa\
    -V gendb:///lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/4_variant_calling/2_GenomicsDB_data_store/\
    -L /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/intervals.list\
    -O BAR_joint_genotyping.vcf

echo "End Joint Genotyping"
#################################################################################################
echo "Printing Statistics"

### COUNTVARIANTS ###
# counts variant records in VCF file regardless of filter status
gatk  --java-options "-Xmx4g" CountVariants\
    -V BAR_joint_genotyping.vcf\
    --output BAR_CountVariants

### VARIANTSTOTABLE ###
#extracts fields from VCF to tab-delimited table
gatk  --java-options "-Xmx4g" VariantsToTable\
    -V BAR_joint_genotyping.vcf\
    --output BAR_VCF.table

echo "Finished Printing Statistics"

module purge
echo "End Job"
