#!/bin/bash

# script to extract reads from CRAM files corresponding to information from VCF file
# shorten prompt
export PROMPT_DIRTRIM=1
# remove headers from VCF file and extract column 3 containing gain/loss ID information
grep -v "##" <SV_vcf_file.vcf>|cut -f 3 > col3_no_header_SV_vcf_file.vcf
# remove column header and replace : and - separators with tab
sed '1d;s/:/\t/g;s-/\t/g' col3_no_header_SV_vcf_file.vcf > col3_tabbed.vcf
# create a bed file with chromosome, start, stop and gain/loss information
awk '{print $3,t,$4,t,$5,t,$2}' col3_tabbed.vcf > SV_vcf.bed
rm col3_no_header_SV_vcf_file.vcf col3_tabbed.vcf
# run samtools to extract reads from CRAM file corresponding to positions in VCF file
samtools view -T GRCh38_full_analysis_set_plus_decoy_hla.fa -L SV_vcf.bed -b -o extracted_positions_cram_vcf.bam cancer_sample_cram_file.cram

