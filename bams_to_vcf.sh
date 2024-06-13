#!/bin/bash

genomefa="/path/to/genome"

### MOVE SCAFFOLD BY SCAFFOLD
bamdir="/path/to/bams/"
vcfdir="/path/to/vcfs/"
gvcfdir="/path/to/gvcfs/"

### CREATE NEEDED DIRECTORIES
###  -p flag will create the full path and will do nothing
###     but exit without error if the path already exists

mkdir -p $vcfdir
mkdir -p $gvcfdir

### CREATE GVCFS FOR EACH SAMPLE

ls $bamdir | grep ".bam" | grep -v "bai" | sed 's/.bam//' | while read bam
do seq 8 | parallel -j 8 \
   java -jar GenomeAnalysisTK.jar -R ${genomefa} \
	-T HaplotypeCaller -I ${bamdir}${bam}.bam \
	--emitRefConfidence GVCF \
	-L CE10_chr{1} \
	-o ${gvcfdir}${bam}.scaffold_{1}.gvcf \
	-variant_index_type LINEAR -variant_index_parameter 128000
done

### MAKE A LIST OF THE CHROMOSOMES YOU WANT TO CREATE VCFS FOR
lgs="1
2
3
4
5
6
7
8
"

### MAKE VCF FOR EACH CHROMOSOME, INCLUDING ALL SITES ON THE CHROMOSOME
for lg in $lgs
do java -jar GenomeAnalysisTK.jar -R $genomefa \
	-T GenotypeGVCFs \
	--includeNonVariantSites \
	$(ls $gvcfdir/*Genome_name${lg}.gvcf | while read v ; do echo -n " --variant $v "; done) \
	-o ${vcfdir}/${lg}.allsites.vcf
   bgzip ${vcfdir}/${lg}.allsites.vcf
   tabix -p vcf VCF/${lg}.allsites.vcf.gz
done

### FILTER VCF TO SITES WITH MINIMAL MISSING DATA
### ('vcftools --max-missing 0.9' allows 10% missing data)

### FILTER TO GENIC REGIONS USING A BED FILE (vcftools --bed))

### USE VCFTOOLS --min-alleles 2 --max-alleles 2 TO FILTER TO BIALLELIC SITES

### PHASE GENOTYPES

beagle="/home/thom_nelson/opt/beagle.27Jan18.7e1.jar"
maxmem="50" # in GB
outpath="/home/colette_berg/sequenceData3.2020/3762/raw/vcfs/phased/"

java -Xss5m -Xmx${maxmem}g -jar $beagle \
     gt=VCF/vcfName.vcf.gz \
     out=${outpath}vcfName.phased \
     nthreads=16 ne=100000
tabix -p vcf ${outpath}vcfName.phased.vcf.gz
