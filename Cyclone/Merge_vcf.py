"""
filename: Merge_vcf.py
description: This script is to integrate the results of genotype filling 
author: HeLei
email: helei1@genomics.cn
last modified: 2023-10-23
"""

import argparse
import os

parser = argparse.ArgumentParser(description='Merge two VCF files.')
parser.add_argument('--vcf_qual_fn', required=True, help='VCF file with QUAL values')
parser.add_argument('--vcf_gp_fn', required=True, help='path of VCF files with GP values')
parser.add_argument('--output', required=True, help='Output VCF file')
parser.add_argument('--bcftools', default= "bcftools", type=str, help="Path of bcftools, bcftools version >= 1.13 is required.")
parser.add_argument('--tabix', default= "tabix", type=str, help="Path of tabix.")
args = parser.parse_args()

def line_iter(vcf):
    with open(vcf,"r") as vr:
        for line in vr:
            if line[0] != "#":
                yield line
            else:
                continue

BCFTOOLS = args.bcftools
TABIX = args.tabix
vcf_qual = args.vcf_qual_fn
vcf_gp = args.vcf_gp_fn

BCFTOOLS = "bcftools"
vcf_qual = "/data2/helei/work/SNP_calling/HG002_gufen/HG002_gufen_chr20_30X/merge_output.vcf.gz"
vcf_gp = "/data2/helei/work/SNP_calling/HG002_gufen/HG002_gufen_chr20_30X/Genotype_filling_final.vcf.gz"
if ".gz" in vcf_qual:
    os.system("{} view --threads 4 -Ov -o {} {}".format(BCFTOOLS,vcf_qual.replace(".gz",""),vcf_qual))
    vcf_qual = vcf_qual.replace(".gz","")
if ".gz" in vcf_gp:
    os.system("{} view --threads 4 -Ov -o {} {}".format(BCFTOOLS,vcf_gp.replace(".gz",""),vcf_gp))
    vcf_gp = vcf_gp.replace(".gz","")

new_header_filter = []
new_header_info = []
new_header_format = []
new_header_contig = []
new_header_line = []
new_header = []
new_header.append('##fileformat=VCFv4.2\n')
new_header.append("##source=LRAPmut v1.0\n")

with open(vcf_gp,"r") as vcf_gp_header:
    for header_line in vcf_gp_header:
        if header_line.startswith('##FILTER='):
            new_header_filter.append(header_line)
        elif header_line.startswith('##INFO='):
            new_header_info.append(header_line)
        elif header_line.startswith('##FORMAT='):
            new_header_format.append(header_line)
        elif header_line[0] != "#":
            break

with open(vcf_qual,"r") as vcf_qual_header:
    for header_line in vcf_qual_header:
        if header_line.startswith('##FILTER='):
            new_header_filter.append(header_line)
        elif header_line.startswith('##INFO='):
            new_header_info.append(header_line)
        elif header_line.startswith('##FORMAT='):
            new_header_format.append(header_line)
        elif header_line.startswith('##contig='):
            new_header_contig.append(header_line)
        elif header_line.startswith('#CHROM'):
            new_header_line.append(header_line)
        elif header_line[0] != "#":
            break

new_header += new_header_filter
new_header += new_header_info
new_header += new_header_format
new_header += new_header_contig
new_header += new_header_line
# Initialize record iterator
iter_qual = line_iter(vcf_qual)
iter_gp = line_iter(vcf_gp)
# Obtain initial records
rec_qual = next(iter_qual, None)
rec_gp = next(iter_gp, None)
with open(args.output, "w") as vcf_out:
    vcf_out.write("".join(new_header))
    while rec_qual is not None and rec_gp is not None:
        if "," in rec_qual.split(":")[-1]:
            rec_qual_af = max([float(i) for i in rec_qual.strip().split(":")[-1].split(",")])
        else:
            rec_qual_af = float(rec_qual.strip().split(":")[-1])
        if (rec_qual_af < 0.1) or (rec_qual.split("\t")[6] != "." and rec_qual.split("\t")[6] != "PASS"):
            rec_qual = next(iter_qual, None)
            continue
        rec_qual_chrom = rec_qual.split("\t")[0]
        if ("_" in rec_qual_chrom) or ("M" in rec_qual_chrom):
            rec_qual = next(iter_qual, None)
            continue
        if "chr" in rec_qual_chrom:
            qual_chrom = rec_qual_chrom.replace("chr","")
        else:
            qual_chrom = rec_qual_chrom
        if qual_chrom not in ["X", "Y"]:
            qual_chrom = int(qual_chrom)
        elif qual_chrom == "X":
            qual_chrom = 23
        elif qual_chrom == "Y":
            qual_chrom = 24
        qual_pos = int(rec_qual.split("\t")[1])
        qual_alts = rec_qual.split("\t")[4]
        qual_ref = rec_qual.split("\t")[3]
        
        rec_gp_af = rec_gp.split("AF=")[1].split(";")[0]
        if "," in rec_gp_af:
            rec_gp_af = max([float(i) for i in rec_gp_af.split(",")])
        else:
            rec_gp_af = float(rec_gp_af)
        if (rec_gp_af < 0.1) or (rec_gp.split("\t")[-1].split(":")[0] == "0/0"):
            rec_gp = next(iter_gp, None)
            continue
        rec_gp_chrom = rec_gp.split("\t")[0]
        if ("_" in rec_gp_chrom) or ("M" in rec_gp_chrom):
            rec_gp = next(iter_gp, None)
            continue
        if "chr" in rec_gp_chrom:
            gp_chrom = rec_gp_chrom.replace("chr","")
        else:
            gp_chrom = rec_gp_chrom
        if gp_chrom not in ["X", "Y"]:
            gp_chrom = int(gp_chrom)
        elif gp_chrom == "X":
            gp_chrom = 23
        elif gp_chrom == "Y":
            gp_chrom = 24
        gp_pos = int(rec_gp.split("\t")[1])
        gp_alts = rec_gp.split("\t")[4]
        gp_ref = rec_gp.split("\t")[3]
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
        if "," in rec_qual.split(":")[-1]:
            rec_qual_af = max([float(i) for i in rec_qual.strip().split(":")[-1].split(",")])
        else:
            rec_qual_af = float(rec_qual.strip().split(":")[-1])
        if (rec_qual_af < 0.1) or (rec_qual.split("\t")[6] != "." and rec_qual.split("\t")[6] != "PASS"):
            rec_qual = next(iter_qual, None)
            continue
        rec_qual_chrom = rec_qual.split("\t")[0]
        if ("_" in rec_qual_chrom) or ("M" in rec_qual_chrom):
            rec_qual = next(iter_qual, None)
            continue
        qual_alts = rec_qual.split("\t")[4]
        qual_ref = rec_qual.split("\t")[3]
        if len(qual_alts) == len(qual_ref):
            vcf_out.write(str(rec_qual))
        rec_qual = next(iter_qual, None)

    while rec_gp is not None:
        rec_gp_af = rec_gp.split("AF=")[1].split(";")[0]
        if "," in rec_gp_af:
            rec_gp_af = max([float(i) for i in rec_gp_af.split(",")])
        else:
            rec_gp_af = float(rec_gp_af)
        if (rec_gp_af < 0.1) or (rec_gp.split("\t")[-1].split(":")[0] == "0/0"):
            rec_gp = next(iter_gp, None)
            continue
        rec_gp_chrom = rec_gp.split("\t")[0]
        if ("_" in rec_gp_chrom) or ("M" in rec_gp_chrom):
            rec_gp = next(iter_gp, None)
            continue
        gp_alts = rec_gp.split("\t")[4]
        gp_ref = rec_gp.split("\t")[3]
        if len(gp_alts) != len(gp_ref):
            vcf_out.write(str(rec_gp))
        rec_gp = next(iter_gp, None)

tmpdir = os.path.join("/".join(args.output.split("/")[0:-1]),"temp")
if not os.path.exists(tmpdir):
    os.system("mkdir -p {}".format(tmpdir))
if vcf_qual.split(".")[-1] == "vcf":
    os.system("rm {}".format(vcf_qual))
if vcf_gp.split(".")[-1] == "vcf":
    os.system("rm {}".format(vcf_gp))
os.system("{} sort -Oz -o {} {} -T {} ; {} -p vcf {} ; rm -r {}".format(BCFTOOLS,args.output+".gz",args.output,tmpdir,TABIX,args.output+".gz",tmpdir))