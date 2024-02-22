"""
filename: bcftools_mpileup.py
description: This script is the preprocessing of genotype filling 
author: HeLei
email: helei1@genomics.cn
last modified: 2023-10-23
"""

import os
import argparse


def bcftools(REF,VCF,CHROM,BAM,TSV,OUT,BCFTOOLS,TABIX):
    os.system("{} mpileup --threads 2 -f {} -I -E -a 'FORMAT/DP' -T {} -r {} {} -Ou | {} call --threads 2 -Aim -C alleles -T {} -Oz -o {}".format(BCFTOOLS,REF, VCF, CHROM, BAM, BCFTOOLS, TSV, OUT))
    os.system("{} -p vcf {}".format(TABIX,OUT))

def vcf_filter(BCFTOOLS,input_vcf,out_putvcf,tabix):
    os.system("{} view --threads 2 -i 'QUAL>0 && TYPE=\"snp\"' -Oz -o {} {}".format(BCFTOOLS,out_putvcf,input_vcf))
    os.system("rm {}; {} {}".format(input_vcf,tabix,out_putvcf))

#-----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Select candidate snp for haplotype filling')
    parser.add_argument('--bam', required=True, type=str, 
                        help='sorted bam file')
    parser.add_argument('--ref', required=True, type=str, 
                        help='reference file')
    parser.add_argument('--bcftools', default= "bcftools", type=str, 
                        help="Path of bcftools, bcftools version >= 1.13 is required.")
    parser.add_argument('--tabix', default= "tabix", type=str, 
                        help="Path of tabix.")
    parser.add_argument('--outdir', required=True, type=str, 
                        help='Output dir')
    parser.add_argument('--chrome', required=True, type=str, 
                        help='Target chrome')
    parser.add_argument('--ref_panel_site_vcf', required=True, type=str, 
                        help='the dir of ref panel vcfs')
    args = parser.parse_args()

    BAM = args.bam
    REF = args.ref
    BCFTOOLS = args.bcftools
    TABIX = args.tabix
    CHROME = args.chrome
    OUT = os.path.join(args.outdir,"{}_bcftools.vcf.gz".format(CHROME))
    REF_PANEL_SITE_VCF = args.ref_panel_site_vcf
    TSV = REF_PANEL_SITE_VCF.replace("sites.vcf.gz","sites.tsv.gz")
    bcftools(REF,REF_PANEL_SITE_VCF,CHROME,BAM,TSV,OUT,BCFTOOLS,TABIX)
    OUT2 = os.path.join(args.outdir,"{}_bcftools_filtered.vcf.gz".format(CHROME))
    vcf_filter(BCFTOOLS,OUT,OUT2,TABIX)

if __name__ == "__main__":
    main()