"""
filename: Genotype_filling_vcf.py
description: This script is the preprocessing of genotype filling 
author: HeLei
email: helei1@genomics.cn
last modified: 2023-10-23
"""

import pysam
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Select candidate snp for haplotype filling')
parser.add_argument('--var_pct_fill_snp', type=float, default=0.1, 
                    help='Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for haplotype filling, default: 0.1')
parser.add_argument('--var_pct_fill_indel', type=float, default=0.2, 
                    help='Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for haplotype filling, default: 0.2')
parser.add_argument('--vcf_fn', required=True,
                    help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file, default: %(default)s")
parser.add_argument('--output_fn', type=str, required=True,
                    help="Define the output folder, required")
parser.add_argument('--chrome', type=str,
                    help="The name of the sequence to be processed")
parser.add_argument('--bcftools', default= "bcftools", type=str,
                    help="Path of bcftools, bcftools version >= 1.13 is required.")
parser.add_argument('--tabix', default= "tabix", type=str,
                    help="Path of tabix.")
args = parser.parse_args()
chrome = args.chrome
candidate_vcf = {}

qual_list_snp = []
qual_list_indel = []
vcffile = args.vcf_fn
BCFTOOLS = args.bcftools
TABIX = args.tabix
with pysam.VariantFile(vcffile,"r") as vcf:
    for rec in vcf:
        CHR = rec.chrom
        if CHR == chrome:
            genotype = str(rec).split('\t')[-1].split(':')[0]
            if genotype == "0/1" or genotype == "1/1":
                if len(rec.ref) == len(rec.alts[0]):
                    qual = rec.qual
                    qual_list_snp.append(qual)
                else:
                    qual = rec.qual
                    qual_list_indel.append(qual)


qual_list_snp = sorted(qual_list_snp)
low_qual_hete_var_pct_snp = 1 - args.var_pct_fill_snp
low_phase_qual_list_snp = qual_list_snp[:int(low_qual_hete_var_pct_snp * len(qual_list_snp))]
qual_cut_off_snp = low_phase_qual_list_snp[-1]

qual_list_indel = sorted(qual_list_indel)
low_qual_hete_var_pct_indel = 1 - args.var_pct_fill_indel
low_phase_qual_list_indel = qual_list_indel[:int(low_qual_hete_var_pct_indel * len(qual_list_indel))]
qual_cut_off_indel = low_phase_qual_list_indel[-1]


outdir = args.output_fn
pileup_vcf_gf = os.path.join(outdir,"{}_pileup_vcf_gf.vcf".format(chrome))
with pysam.VariantFile(args.vcf_fn,"r") as invcf, pysam.VariantFile(pileup_vcf_gf,"w",header=invcf.header) as outvcf:
     for rec in invcf:
        CHR = rec.chrom
        if CHR == chrome:
            genotype = str(rec).split('\t')[-1].split(':')[0]
            if genotype == "0/1" or genotype == "1/1":
                if len(rec.ref) == len(rec.alts[0]):
                    qual = rec.qual
                    if qual >= qual_cut_off_snp:
                        outvcf.write(rec)
                else:
                    qual = rec.qual
                    if qual >= qual_cut_off_indel:
                        outvcf.write(rec)

subprocess.run("{} sort {} -Oz -o {} -T {} ; {} -p vcf {} ; rm {}".format(BCFTOOLS, pileup_vcf_gf,pileup_vcf_gf+".gz", outdir, TABIX, pileup_vcf_gf+".gz",pileup_vcf_gf),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)