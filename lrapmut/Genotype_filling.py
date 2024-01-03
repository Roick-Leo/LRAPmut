"""
filename: Genotype_filling.py
description: This is the script for genotype filling
author: HeLei
email: helei1@genomics.cn
last modified: 2023-10-23
"""

import os
import argparse
import subprocess
import sys

def glimpse(GLIMPSE_phase_static, GLIMPSE_ligate_static, outdir, chrome, chunks, vcffile, refpanel, refmap, BCFTOOLS, TABIX):
    listfile = "{}/{}.txt".format(outdir+"/tmp", chrome + '_pileup_list')
    with open(chunks, "r") as chunkline, open(listfile, "w") as lf:
        for line in chunkline:
            ID = line.split("\t")[0]
            irg = line.split("\t")[2]
            org = line.split("\t")[3]
            out = "{}/{}.bcf".format(outdir+"/tmp",chrome + '_' + ID)
            commad_GLIMPSE_phase_static = [GLIMPSE_phase_static,
                                           "--thread", "2",
                                           "--input", vcffile,
                                           "--reference", refpanel,
                                           "--map", refmap,
                                           "--impute-reference-only-variants",
                                           "--input-region", irg,
                                           "--output-region", org,
                                           "--output",out]
            try:
                subprocess.check_call(commad_GLIMPSE_phase_static)
                subprocess.run("{} index -f {}".format(BCFTOOLS,out),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                lf.write(out + "\n")
            except subprocess.CalledProcessError as e:
                print("[INFO]Error in genotype_filling step1:", e.returncode)
                sys.exit()

    merge_bcf = "{}/{}.bcf".format(outdir, chrome + '_merged')
    commad_GLIMPSE_ligated = [GLIMPSE_ligate_static,
                              "--thread", "2",
                              "--input", listfile,
                              "--output", merge_bcf,
                              ]
    try:
        subprocess.check_call(commad_GLIMPSE_ligated)
    except subprocess.CalledProcessError as e:
        print("[INFO]Error in genotype_filling step2:", e.returncode)
        sys.exit()
    subprocess.run("{} index --threads 2 -f {}".format(BCFTOOLS,merge_bcf),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    finalvcf = "{}/{}.vcf.gz".format(outdir, chrome)
    subprocess.run("{} view --thread 2 {} -Oz -o {} ; {} -p vcf {} ; rm {}*".format(BCFTOOLS, merge_bcf, finalvcf, TABIX, finalvcf, merge_bcf),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)

#########################################################################################################################################################
def main():
    parser = argparse.ArgumentParser(description='Bcftools-glimpse pipline')
    parser.add_argument('--vcf_fn', required=True,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file, default: %(default)s")
    parser.add_argument('--ref_panel_fn', required=True,
                        help='the dir of ref panel vcfs')
    parser.add_argument('--reference_map', required=True,
                        help='the file of reference maps')
    parser.add_argument('--chunks', required=True,
                        help='path of chunks')
    parser.add_argument('--outdir', required=True,
                        help='Output dir')
    parser.add_argument('--chrome', required=True,
                        help='out put prefix')
    parser.add_argument('--glimpse', required=True,
                        help='the dir of glimpse')
    parser.add_argument('--bcftools', default= "bcftools", type=str, 
                        help="Path of bcftools, bcftools version >= 1.13 is required.")
    parser.add_argument('--tabix', default= "tabix", type=str, 
                        help="Path of tabix.")
    args = parser.parse_args()

    glimpse_dir = args.glimpse
    GLIMPSE_phase_static = os.path.join(glimpse_dir,"GLIMPSE_phase_static")
    GLIMPSE_ligate_static = os.path.join(glimpse_dir,"GLIMPSE_ligate_static")
    refdir = args.ref_panel_fn
    reference_map_dir = args.reference_map
    chrome = args.chrome
    outdir = args.outdir
    chunks = args.chunks
    vcffile = args.vcf_fn
    BCFTOOLS = args.bcftools
    TABIX = args.tabix

    glimpse(GLIMPSE_phase_static, GLIMPSE_ligate_static, outdir, chrome, chunks, vcffile, refdir, reference_map_dir, BCFTOOLS, TABIX)

if __name__ == "__main__":
    main()