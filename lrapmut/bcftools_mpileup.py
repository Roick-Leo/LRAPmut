"""
filename: bcftools_mpileup.py
description: This script is the preprocessing of genotype filling 
author: HeLei
email: helei1@genomics.cn
last modified: 2023-10-23
"""

import os
import pysam
import argparse

def bcftools(REF,VCF,CHROM,BAM,TSV,OUT,BCFTOOLS,TABIX):
    os.system("{} mpileup --threads 2 -f {} -I -E -a 'FORMAT/DP' -T {} -r {} {} -Ou | {} call --threads 2 -Aim -C alleles -T {} -Oz -o {}".format(BCFTOOLS,REF, VCF, CHROM, BAM, BCFTOOLS, TSV, OUT))
    os.system("{} -p vcf {}".format(TABIX,OUT))

def vcf_filter(input_vcf,out_putvcf,tabix):
    with pysam.VariantFile(input_vcf,"r") as vcfin, pysam.VariantFile(out_putvcf,"w",header=vcfin.header) as vcfout:
        for rec in vcfin:
            qual = rec.qual
            if isinstance(qual,float):
                if rec.qual > 0 and len(rec.ref) == len(rec.alts[0]):
                    vcfout.write(rec)
    os.system("bcftools view --threads 2 -Oz -o {} {}".format(out_putvcf+".gz",out_putvcf))
    os.system("rm {} {}; {} {}".format(input_vcf,out_putvcf,tabix,out_putvcf+".gz"))

#-----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Select candidate snp for haplotype filling')
    parser.add_argument('--bam', required=True, 
                        help='sorted bam file')
    parser.add_argument('--ref', required=True, 
                        help='reference file')
    parser.add_argument('--bcftools', default= "bcftools", type=str, 
                        help="Path of bcftools, bcftools version >= 1.13 is required.")
    parser.add_argument('--tabix', default= "tabix", type=str, 
                        help="Path of tabix.")
    parser.add_argument('--outdir', required=True,
                        help='Output dir')
    parser.add_argument('--chrome', required=True,
                        help='Target chrome')
    parser.add_argument('--ref_panel_site_vcf', required=True,
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
    OUT2 = os.path.join(args.outdir,"{}_bcftools_filtered.vcf".format(CHROME))
    vcf_filter(OUT,OUT2,TABIX)

if __name__ == "__main__":
    main()