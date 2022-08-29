#!/bin/bash

module load bioinfo/vcftools/0.1.16

mydir=$PWD
scratch=$(mktemp -dp /scratch)

vcf_file=vietcaf_general_DP_missing_indv_singletons_Filtered_biallelic.recode.vcf
vcf_out=vietcaf_final_chr01_C013_C034
cd $scratch

# transfer input file to scratch
rsync -vauP $mydir/$vcf_file ./

# extract snps by chromosome
vcftools \
--vcf $vcf_file \
--chr CC1.8.Chr01 \
--indv "C013-MERGED" --indv "C034-MERGED" \
--out $vcf_out \
--recode


# transfer output files to nas
rsync -vauP $vcf_out.recode.vcf $mydir

cd ../
rm -rf $scratch