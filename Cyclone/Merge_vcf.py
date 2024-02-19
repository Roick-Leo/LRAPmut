"""
filename: Merge_vcf.py
description: This script is to integrate the results of genotype filling 
author: HeLei
email: helei1@genomics.cn
last modified: 2023-10-23
"""

import pysam
import argparse
import os

parser = argparse.ArgumentParser(description='Merge two VCF files.')
parser.add_argument('--vcf_qual_fn', required=True, help='VCF file with QUAL values')
parser.add_argument('--vcf_gp_fn', required=True, help='path of VCF files with GP values')
parser.add_argument('--output', required=True, help='Output VCF file')
parser.add_argument('--bcftools', default= "bcftools", type=str, help="Path of bcftools, bcftools version >= 1.13 is required.")
parser.add_argument('--tabix', default= "tabix", type=str, help="Path of tabix.")
args = parser.parse_args()

# Loading vcf file
BCFTOOLS = args.bcftools
TABIX = args.tabix
vcf_qual = pysam.VariantFile(args.vcf_qual_fn,"r")
vcf_gp = pysam.VariantFile(args.vcf_gp_fn,"r")
new_header_filter = pysam.VariantHeader()
new_header_info = pysam.VariantHeader()
new_header_format = pysam.VariantHeader()
new_header_contig = pysam.VariantHeader()
new_header = pysam.VariantHeader()
new_header.add_line('##fileformat=VCFv4.2')
new_header.add_line("##source=LRAPmut v0.1")
new_header.add_sample("SAMPLE")
for record in vcf_gp.header.records:
    if (record.key == 'FILTER'):
        new_header_filter.add_record(record)
    elif (record.key == 'INFO'):
        new_header_info.add_record(record)
    elif record.key == 'FORMAT':
        new_header_format.add_record(record)
    elif record.key == 'contig':
        new_header_contig.add_record(record)

for record in vcf_qual.header.records:
    if (record.key == 'FILTER'):
        new_header_filter.add_record(record)
    elif (record.key == 'INFO'):
        new_header_info.add_record(record)
    elif record.key == 'FORMAT':
        new_header_format.add_record(record)
    elif record.key == 'contig':
        new_header_contig.add_record(record)

for record in new_header_filter.records:
    new_header.add_record(record)
for record in new_header_info.records:
    new_header.add_record(record)
for record in new_header_format.records:
    new_header.add_record(record)
for record in new_header_contig.records:
    new_header.add_record(record)

vcf_out = pysam.VariantFile(args.output, "w", header=new_header)
vcf_out.close()
# Initialize record iterator
iter_qual = iter(vcf_qual)
iter_gp = iter(vcf_gp)

# Obtain initial records
rec_qual = next(iter_qual, None)
rec_gp = next(iter_gp, None)
with open(args.output, "a") as vcf_out:
    while rec_qual is not None and rec_gp is not None:
        if (rec_gp.info["AF"][0] < 0.1) or (str(rec_gp).split("\t")[-1].split(":")[0] == "0/0"):
            rec_gp = next(iter_gp, None)
            continue
        if not isinstance(rec_qual.samples[rec_qual.samples.keys()[0]]["AF"],float):
            rec_qual_af = max(rec_qual.samples[rec_qual.samples.keys()[0]]["AF"])
        else:
            rec_qual_af = rec_qual.samples[rec_qual.samples.keys()[0]]["AF"]
        if (rec_qual_af < 0.1) or (rec_qual.filter.keys()[0]!= "." and rec_qual.filter.keys()[0]!= "PASS"):
            rec_qual = next(iter_qual, None)
            continue
        if ("_" in rec_qual.chrom) or ("chrM" in rec_qual.chrom):
            rec_qual = next(iter_qual, None)
            continue
        if ("_" in rec_gp.chrom) or ("chrM" in rec_gp.chrom):
            rec_gp = next(iter_gp, None)
            continue
        qual_chrom = rec_qual.chrom.replace("chr","")
        if qual_chrom not in ["X", "Y"]:
            qual_chrom = int(qual_chrom)
        elif qual_chrom == "X":
            qual_chrom = 23
        elif qual_chrom == "Y":
            qual_chrom = 24
        qual_pos = int(rec_qual.pos)
        qual_alts = rec_qual.alts[0]
        qual_ref = rec_qual.ref
        gp_chrom = rec_gp.chrom.replace("chr","")
        if gp_chrom not in ["X", "Y"]:
            gp_chrom = int(gp_chrom)
        elif gp_chrom == "X":
            gp_chrom = 23
        elif gp_chrom == "Y":
            gp_chrom = 24
        gp_pos = int(rec_gp.pos)
        gp_alts = rec_gp.alts[0]
        gp_ref = rec_gp.ref
        if (qual_chrom, qual_pos) == (gp_chrom, gp_pos):
            if qual_alts == gp_alts:
                if len(qual_alts) != len(qual_ref):
                    vcf_out.write(str(rec_gp))
                else:
                    vcf_out.write(str(rec_qual))
            elif len(gp_alts) != len(gp_ref):
                vcf_out.write(str(rec_gp))
            elif len(qual_alts) == len(qual_ref):
                vcf_out.write(str(rec_qual))
            rec_qual = next(iter_qual, None)
            rec_gp = next(iter_gp, None)
        elif (qual_chrom, qual_pos) < (gp_chrom, gp_pos):
            if len(qual_alts) == len(qual_ref):
                vcf_out.write(str(rec_qual))
            rec_qual = next(iter_qual, None)
        else:
            if len(gp_alts) != len(gp_ref):
                vcf_out.write(str(rec_gp))
            rec_gp = next(iter_gp, None)

    while rec_qual is not None:
        if not isinstance(rec_qual.samples[rec_qual.samples.keys()[0]]["AF"],float):
            rec_qual_af = max(rec_qual.samples[rec_qual.samples.keys()[0]]["AF"])
        else:
            rec_qual_af = rec_qual.samples[rec_qual.samples.keys()[0]]["AF"]
        if (rec_qual_af < 0.1) or (rec_qual.filter.keys()[0]!= "." and rec_qual.filter.keys()[0]!= "PASS"):
            rec_qual = next(iter_qual, None)
            continue
        qual_alts = rec_qual.alts[0]
        qual_ref = rec_qual.ref
        if len(qual_alts) == len(qual_ref):
            vcf_out.write(str(rec_qual))
        rec_qual = next(iter_qual, None)

    while rec_gp is not None:
        if (rec_gp.info["AF"][0] < 0.1) or (str(rec_gp).split("\t")[-1].split(":")[0] == "0/0"):
            rec_gp = next(iter_gp, None)
            continue
        gp_ref = rec_gp.ref
        gp_alts = rec_gp.alts[0]
        if len(gp_alts) != len(gp_ref):
            vcf_out.write(str(rec_gp))
        rec_gp = next(iter_gp, None)

vcf_qual.close()
vcf_gp.close()
tmpdir = os.path.join("/".join(args.output.split("/")[0:-1]),"temp")
if not os.path.exists(tmpdir):
    os.system("mkdir -p {}".format(tmpdir))
os.system("{} sort -Oz -o {} {} -T {} ; {} -p vcf {} ; rm -r {}".format(BCFTOOLS,args.output+".gz",args.output,tmpdir,TABIX,args.output+".gz",tmpdir))