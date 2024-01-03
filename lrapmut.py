"""
filename: Merge_vcf.py
description: This is the main script of long-reads snp calling
author: HeLei
email: helei1@genomics.cn
last modified: 2023-10-23
"""

import os
import argparse
import subprocess
import pysam
import sys
from datetime import datetime

# Get script location
script_path = os.path.dirname(os.path.abspath(__file__))

# The options
parser = argparse.ArgumentParser(description='whole calling workflow of cyclone')
parser.add_argument('-b', '--bam_fn', required=True, type=str, help="BAM file input. The input file must be samtools indexed.")
parser.add_argument('-f', '--ref_fn', required=True, type=str, help="FASTA reference file input. The input file must be samtools indexed.")
parser.add_argument('-t', '--threads', default=4, required=True, type=int, help="Max #threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The #chunks being processed simultaneously is ceil(#threads/4)*3. 3 is the overloading factor.")
parser.add_argument('-m', '--model_path', required=True, type=str, help="The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00002, pileup.data-00001-of-00002 pileup.index, full_alignment.data-00000-of-00002, full_alignment.data-00001-of-00002 and full_alignment.index).")
parser.add_argument('-p', '--platform', required=True, type=str, default="ont", help="Select the sequencing platform of the input. Possible options: {ont,hifi,ilmn}.")
parser.add_argument('-o', '--output', required=True, type=str, help="VCF/GVCF output directory.")
parser.add_argument('-refpanel', '--reference_panel_fn', required=True, type=str, help="Path of reference panel.")
parser.add_argument('-refmap', '--reference_map_fn', required=True, type=str, help="Path of reference maps.")
parser.add_argument('--glimpse_fn', required=True, type=str, help="Path of glimpse.")
parser.add_argument('--bed_fn', default="EMPTY", type=str, help="Call variants only in the provided bed regions.")
parser.add_argument('--vcf_fn', default="EMPTY", type=str, help="Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.")
parser.add_argument('--ctg_name', default="EMPTY", help="The name of the sequence to be processed.")
parser.add_argument('--sample_name', default= "SAMPLE", type=str, help="Define the sample name to be shown in the VCF file.")
parser.add_argument('--chunk_num', default=0, type=int, help="CHUNK_NUM")
parser.add_argument('--chunk_size', default=5000000, type=int, help="The size of each chuck for parallel processing, default: 5000000.")
parser.add_argument('--qual', default=2, type=int, help="If set, variants with >$qual will be marked PASS, or LowQual otherwise.")
parser.add_argument('--samtools', default= "samtools", type=str, help="Path of samtools, samtools version >= 1.10 is required.")
parser.add_argument('--bcftools', default= "bcftools", type=str, help="Path of bcftools, bcftools version >= 1.13 is required.")
parser.add_argument('--python', default= "python3", type=str, help="Path of python.")
parser.add_argument('--tabix', default= "tabix", type=str, help="Path of tabix.")
parser.add_argument('--pypy', default= "pypy3", type=str, help="Path of pypy3, pypy3 >= 3.6 is required.")
parser.add_argument('--parallel', default= 'parallel', type=str, help="Path of parallel, parallel >= 20191122 is required.")
parser.add_argument('--whatshap', default= "whatshap", type=str, help="Path of whatshap, whatshap >= 1.0 is required.")
parser.add_argument('--longphase', default= "longphase", type=str, help="Path of longphase, longphase >= 1.0 is required.")
parser.add_argument('--var_pct_full', default=0.7, type=float, help="EXPERIMENTAL: Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.7.")
parser.add_argument('--ref_pct_full', default=0.1, type=float, help="EXPERIMENTAL: Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for ilmn and hifi, 0.1 for ont.")
parser.add_argument('--var_pct_phasing', default=0.7, type=float, help="EXPERIMENTAL: Specify an expected percentage of high quality 0/1 variants used in WhatsHap phasing, default: 0.8 for ont guppy5 and 0.7 for other platforms.")
parser.add_argument('--snp_min_af', default=0.08, type=float, help="Minimum SNP AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default:0.08.")
parser.add_argument('--indel_min_af', default=0.15, type=float, help="Minimum Indel AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default:0.15")
parser.add_argument('--min_mq', default=5, type=int, help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: 5.")
parser.add_argument('--min_coverage', default=2, type=int, help="EXPERIMENTAL: Minimum coverage required to call a variant, default: 2.")
parser.add_argument('--min_contig_size', default=0, type=int, help="EXPERIMENTAL: If set, contigs with contig size<$min_contig_size are filtered, default: 0.")
parser.add_argument('--pileup_model_prefix', default="pileup", type=str, help="EXPERIMENTAL: Model prefix in pileup calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index. default: pileup.")
parser.add_argument('--fa_model_prefix', default="full_alignment", type=str, help="EXPERIMENTAL: Model prefix in full-alignment calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index, default: full_alignment.")
parser.add_argument('--gvcf', action='store_true', help="Enable GVCF output, default: disable.")
parser.add_argument('--pileup_only', action='store_true', help="Use the pileup model only when calling, default: disable.")
parser.add_argument('--fast_mode', action='store_true', help="EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: disable.")
parser.add_argument('--call_snp_only', action='store_true', help="EXPERIMENTAL: Call candidates pass SNP minimum AF only, ignore Indel candidates, default: disable.")
parser.add_argument('--print_ref_calls', action='store_true', help="Show reference calls (0/0) in VCF file, default: disable.")
parser.add_argument('--haploid_precise', action='store_true', help="EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant, default: disable.")
parser.add_argument('--haploid_sensitive', action='store_true', help="EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant, default: disable.")
parser.add_argument('--include_all_ctgs', action='store_true', help="Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}, default: disable.")
parser.add_argument('--no_phasing_for_fa', action='store_true', help="EXPERIMENTAL: Call variants without whatshap phasing in full alignment calling, default: disable.")
parser.add_argument('--remove_intermediate_dir', action='store_true', help="Remove intermediate directory, including intermediate phased BAM, pileup and full-alignment results. default: disable.")
parser.add_argument('--use_whatshap_for_intermediate_phasing', action='store_true', help="Phase high-quality heterozygous variants using whatshap for full-alignment model calling, default: enable.")
parser.add_argument('--use_longphase_for_intermediate_phasing', action='store_true', help="Phase high-quality heterozygous variants using longphase for full-alignment model calling, default: disable.")
parser.add_argument('--use_whatshap_for_final_output_phasing', action='store_true', help="Phase the output variants using whatshap, default: disable.")
parser.add_argument('--use_longphase_for_final_output_phasing', action='store_true', help="Phase the output variants using longphase, default: disable.")
parser.add_argument('--use_whatshap_for_final_output_haplotagging', action='store_true', help="Haplotag input BAM using output phased variants using whatshap, default: disable.")
parser.add_argument('--enable_phasing', action='store_true', help="It means `--use_whatshap_for_final_output_phasing`. The option is retained for backward compatibility.")
parser.add_argument('--enable_long_indel', action='store_true', help="EXPERIMENTAL: Call long Indel variants(>50 bp), default: disable.")
parser.add_argument('--keep_iupac_bases', action='store_true', help="EXPERIMENTAL: Keep IUPAC reference and alternate bases, default: convert all IUPAC bases to N.")
parser.add_argument('--use_gpu', action='store_true', help="use gpu or not")
args = parser.parse_args()

BCFTOOLS = args.bcftools
TABIX = args.tabix
PYTHON = args.python
CLAIR3 = os.path.join(script_path,"clair3.py")
OUTPUT_FOLDER = args.output
if not os.path.exists(OUTPUT_FOLDER):
    os.system("mkdir -p {}".format(OUTPUT_FOLDER))
MODEL_PATH = args.model_path
PILEUP_PREFIX = args.pileup_model_prefix
FA_PREFIX = args.fa_model_prefix
PILEUP_CHECKPOINT_PATH = os.path.join(MODEL_PATH,PILEUP_PREFIX)
FULL_ALIGNMENT_CHECKPOINT_PATH=os.path.join(MODEL_PATH,FA_PREFIX)
LOG_PATH=os.path.join(OUTPUT_FOLDER,"log")
TMP_FILE_PATH=os.path.join(OUTPUT_FOLDER,"tmp")
SPLIT_BED_PATH=os.path.join(TMP_FILE_PATH,"split_beds")
PILEUP_VCF_PATH=os.path.join(TMP_FILE_PATH,"pileup_output")
GVCF_TMP_PATH=os.path.join(TMP_FILE_PATH,"gvcf_tmp_output")
PHASE_OUTPUT_PATH=os.path.join(TMP_FILE_PATH,"phase_output")
FULL_ALIGNMENT_OUTPUT_PATH=os.path.join(TMP_FILE_PATH,"full_alignment_output")
PHASE_VCF_PATH=os.path.join(PHASE_OUTPUT_PATH,"phase_vcf")
PHASE_BAM_PATH=os.path.join(PHASE_OUTPUT_PATH,"phase_bam")
CANDIDATE_BED_PATH=os.path.join(FULL_ALIGNMENT_OUTPUT_PATH,"candidate_bed")
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['GOTO_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
print("[{}] Check environment variables".format(formatted_time))

BAM_FILE_PATH = args.bam_fn
if args.bed_fn == "EMPTY":
    BED_FILE_PATH=""
else:
    BED_FILE_PATH = args.bed_fn
REFERENCE_FILE_PATH = args.ref_fn
VCF_FILE_PATH = args.vcf_fn
CONTIGS = args.ctg_name
CHUNK_NUM = str(args.chunk_num)
CHUNK_SIZE = str(args.chunk_size)
INCLUDE_ALL_CTGS = str(args.include_all_ctgs)
THREADS = args.threads
PYPY = args.pypy
SAMTOOLS = args.samtools
WHATSHAP = args.whatshap
PARALLEL = args.parallel
QUAL = str(args.qual)
SAMPLE = args.sample_name
PRO = str(args.var_pct_full)
REF_PRO = str(args.ref_pct_full)
SNP_AF = str(args.snp_min_af)
INDEL_AF = str(args.indel_min_af)
MIN_CONTIG_SIZE = str(args.min_contig_size)
CMD_FN = os.path.join(TMP_FILE_PATH,"CMD")

command_CheckEnvs = [
    PYTHON, CLAIR3, "CheckEnvs",
    "--bam_fn", BAM_FILE_PATH,
    "--bed_fn", BED_FILE_PATH,
    "--output_fn_prefix", OUTPUT_FOLDER,
    "--ref_fn", REFERENCE_FILE_PATH,
    "--vcf_fn", VCF_FILE_PATH,
    "--ctg_name", CONTIGS,
    "--chunk_num", CHUNK_NUM,
    "--chunk_size", CHUNK_SIZE,
    "--include_all_ctgs", INCLUDE_ALL_CTGS,
    "--threads", str(THREADS),
    "--python", PYTHON,
    "--pypy", PYPY,
    "--samtools", SAMTOOLS,
    "--whatshap", WHATSHAP,
    "--parallel", PARALLEL,
    "--qual", QUAL,
    "--sampleName", SAMPLE,
    "--var_pct_full", PRO,
    "--ref_pct_full", REF_PRO,
    "--snp_min_af", SNP_AF,
    "--indel_min_af", INDEL_AF,
    "--min_contig_size", MIN_CONTIG_SIZE,
    "--cmd_fn", CMD_FN
]
# try:
#     subprocess.check_call(command_CheckEnvs)
#     print("[INFO]Complete environmental checking")
# except subprocess.CalledProcessError as e:
#     print("[INFO]Error in environmental checking:", e.returncode)
#     sys.exit()

CHR = []
with open("{}/CONTIGS".format(TMP_FILE_PATH),"r") as ct:
    for line in ct:
        CHR.append(line.strip())
if len(CHR) == 0:
    print("[INFO] Exit in environment checking")
    sys.exit()

# pileup using BiLSTM
#-----------------------------------------------------------------------------------------------------------------------
import concurrent.futures
def run_subprocess(command):
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError as e:
        return("[Error]:", e.returncode)
    # subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

chunk_list = os.path.join(TMP_FILE_PATH,"CHUNK_LIST")
THREADS_LOW = int(THREADS)
LONGPHASE_THREADS = int(THREADS)
if THREADS_LOW < 1:
    THREADS_LOW = 1
if LONGPHASE_THREADS < 1:
    LONGPHASE_THREADS = 1
PLATFORM = args.platform
if PLATFORM == "ont":
    LP_PLATFORM = "ont"
else:
    LP_PLATFORM = "pb"

now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
print("[{}] STEP1 Call variants using pileup model".format(formatted_time))
MIN_MQ = str(args.min_mq)
FAST_MODE = str(args.fast_mode)
MIN_COV = str(args.min_coverage)
SNP_ONLY = str(args.call_snp_only)
GVCF = str(args.gvcf)
ENABLE_LONG_INDEL = str(args.enable_long_indel)
KEEP_IUPAC_BASES = str(args.keep_iupac_bases)
USE_GPU = str(args.use_gpu)
with concurrent.futures.ThreadPoolExecutor(THREADS) as executor, open(chunk_list,"r") as cl:
    for line in cl:
        linelist = line.strip().split(" ")
        CALL_FN = os.path.join(PILEUP_VCF_PATH,"pileup_{}_{}.vcf".format(linelist[0],linelist[1]))
        EXTEND_BED = os.path.join(SPLIT_BED_PATH,"{}".format(linelist[0]))
        command = [PYTHON, CLAIR3, 'CallVariantsFromCffi',
            '--chkpnt_fn', PILEUP_CHECKPOINT_PATH,
            '--bam_fn', BAM_FILE_PATH,
            '--call_fn', CALL_FN,
            '--sampleName', SAMPLE,
            '--ref_fn', REFERENCE_FILE_PATH,
            '--extend_bed', EXTEND_BED,
            '--bed_fn', BED_FILE_PATH,
            '--vcf_fn', VCF_FILE_PATH,
            '--ctgName', linelist[0],
            '--chunk_id', linelist[1],
            '--chunk_num', linelist[2],
            '--platform', PLATFORM,
            '--fast_mode', FAST_MODE,
            '--snp_min_af', SNP_AF,
            '--indel_min_af', INDEL_AF,
            '--minMQ', MIN_MQ,
            '--minCoverage', MIN_COV,
            '--call_snp_only', SNP_ONLY,
            '--gvcf', GVCF,
            '--enable_long_indel', ENABLE_LONG_INDEL,
            '--samtools', SAMTOOLS,
            '--temp_file_dir', GVCF_TMP_PATH,
            '--pileup',
            '--keep_iupac_bases', KEEP_IUPAC_BASES,
            '--use_gpu', USE_GPU
        ]
        # executor.submit(run_subprocess,command)

command_pileup_SortVcf = [PYPY, CLAIR3, "SortVcf",
    "--input_dir", PILEUP_VCF_PATH,
    "--vcf_fn_prefix", "pileup",
    "--output_fn", "{}/pileup.vcf".format(OUTPUT_FOLDER),
    "--sampleName", SAMPLE,
    "--ref_fn", REFERENCE_FILE_PATH,
    "--contigs_fn", "{}/CONTIGS".format(TMP_FILE_PATH),
    "--cmd_fn", CMD_FN
]
# subprocess.run(command_pileup_SortVcf)

pileup_vcf = "{}/pileup.vcf.gz".format(OUTPUT_FOLDER)
line_num = 0
with pysam.VariantFile(pileup_vcf,"r") as pv:
    for rec in pv:
        line_num += 1
        if line_num > 100:
            break
if line_num < 2:
    print("[INFO] Exit in pileup variant calling")
    sys.exit()

if args.pileup_only == True:
    if args.remove_intermediate_dir == True:
        print("[INFO] Removing intermediate files in {}/tmp".format(OUTPUT_FOLDER))
        subprocess.run("rm -rf {}/tmp".format(OUTPUT_FOLDER), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("[INFO] Only call pileup output with --pileup_only, output file: {}/pileup.vcf.gz".format(OUTPUT_FOLDER))
    print("[INFO] Finish calling!")
    sys.exit()

# Genotype filling
#-----------------------------------------------------------------------------------------------------------------------
now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
print("[{}] STEP2 Genotype filling".format(formatted_time))
Genotype_filling_fn = os.path.join(TMP_FILE_PATH,"Genotype_filling_output")
Genotype_filling_fn_tmp = os.path.join(TMP_FILE_PATH,"Genotype_filling_output/tmp")
if not os.path.exists(Genotype_filling_fn_tmp):
    subprocess.run("mkdir -p {}".format(Genotype_filling_fn_tmp),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)

refpanel = args.reference_panel_fn
refpanel_dict = {}
for root,dirs,files in os.walk(refpanel):
    for file in files:
        if file.split("_")[-1] == "panel.vcf.gz":
            key = file.split(".")[1]
            file = os.path.join(root,file)
            refpanel_dict[key] = file

refmap = args.reference_map_fn
refmap_dict = {}
for root,dirs,files in os.walk(refmap):
    for file in files:
        if ("b38.gmap.gz" in file) and ("_" not in file):
            key = file.split(".")[0]
            file = os.path.join(root,file)
            refmap_dict[key] = file

Preprocess_genotype_filling_vcf = os.path.join(script_path,"Cyclone/bcftools_mpileup.py")
with concurrent.futures.ThreadPoolExecutor(THREADS) as executor:
    for chrome in CHR:
        if chrome in refpanel_dict.keys():
            REF_PANEL_SITE_VCF = refpanel_dict[chrome].replace("panel.vcf.gz","panel.sites.vcf.gz")
            command_Genotype_filling_vcf = [PYTHON, Preprocess_genotype_filling_vcf,
                                            "--bam", BAM_FILE_PATH,
                                            "--ref", REFERENCE_FILE_PATH,
                                            "--bcftools", BCFTOOLS,
                                            "--tabix", TABIX, 
                                            "--outdir", Genotype_filling_fn_tmp,
                                            "--chrom", chrome,
                                            "--ref_panel_site_vcf",REF_PANEL_SITE_VCF]
            os.system(" ".join(command_Genotype_filling_vcf))
            # executor.submit(run_subprocess,command_Genotype_filling_vcf)

Genotype_filling = os.path.join(script_path,"Cyclone/Genotype_filling.py")
glimpse_dir = args.glimpse_fn
with concurrent.futures.ThreadPoolExecutor(int(THREADS/2)) as executor:
    for chrome in CHR:
        if chrome in refpanel_dict.keys():
            panel = refpanel_dict[chrome]
            refmap_file = refmap_dict[chrome]
            chunk_txt = refpanel+"/chunks."+chrome+".txt"
            vcf_file = os.path.join(Genotype_filling_fn_tmp,"{}_bcftools_filtered.vcf.gz".format(chrome))
            command_Genotype_filling_vcf = [PYTHON, Genotype_filling,
                                            "--vcf_fn", vcf_file,
                                            "--ref_panel_fn", panel,
                                            "--reference_map", refmap_file,
                                            "--outdir", Genotype_filling_fn,
                                            "--chrome", chrome,
                                            "--glimpse",glimpse_dir,
                                            "--chunks",chunk_txt,
                                            "--bcftools",BCFTOOLS,
                                            "--tabix", TABIX]
            executor.submit(run_subprocess,command_Genotype_filling_vcf)

# Whatshap phasing and haplotaging
#-----------------------------------------------------------------------------------------------------------------------
non_phase = args.no_phasing_for_fa
PHASING_PCT = args.var_pct_phasing 
if non_phase == True:
    now = datetime.now()
    formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
    print("[{}] STEP3 No phasing for full alignment calling".format(formatted_time))
    for chrome in CHR:
        subprocess.run('ln -sf {} {}.bam'.format(BAM_FILE_PATH,os.path.join(PHASE_BAM_PATH,CHR)), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if os.path.isfile(os.path.splitext(BAM_FILE_PATH)[0] + ".bai"):
            subprocess.run("ln -sf {}.bai {}/{}.bam.bai".format(os.path.splitext(BAM_FILE_PATH)[0],PHASE_BAM_PATH,CHR), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
else:
    now = datetime.now()
    formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
    print("[{}] STEP3 Select heterozygous SNP variants for Whatshap phasing and haplotagging".format(formatted_time))
    subprocess.run("gzip -fdc {}/pileup.vcf.gz | {} {} SelectQual --phase --output_fn {} --var_pct_phasing {}".format(OUTPUT_FOLDER,PYPY,CLAIR3,PHASE_VCF_PATH,PHASING_PCT), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with concurrent.futures.ThreadPoolExecutor(THREADS) as executor:
        for chrome in CHR:
            command_SelectHetSnp = [PYPY, CLAIR3,"SelectHetSnp",
            "--vcf_fn","{}/pileup.vcf.gz".format(OUTPUT_FOLDER),
            "--split_folder", PHASE_VCF_PATH,
            "--ctgName", chrome]
            executor.submit(run_subprocess,command_SelectHetSnp)
    USE_LONGPHASE = args.use_longphase_for_intermediate_phasing
    if USE_LONGPHASE:
        now = datetime.now()
        formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
        print("[{}] STEP4 Phase VCF file using LongPhase".format(formatted_time))
        LONGPHASE = args.longphase
        with concurrent.futures.ThreadPoolExecutor(THREADS) as executor:
            for chrome in CHR:
                command_longphase = [LONGPHASE, "phase",
                            "-s", "{}/{}.vcf".format(PHASE_VCF_PATH,chrome),
                            "-b", BAM_FILE_PATH,
                            "-r", REFERENCE_FILE_PATH,
                            "-t", LONGPHASE_THREADS,
                            "-o", "{}/phased_{}".format(PHASE_VCF_PATH,chrome),
                            "--{}".format(LP_PLATFORM)]
                executor.submit(run_subprocess,command_longphase)
    else:
        now = datetime.now()
        formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
        print("[{}] STEP4 Phase VCF file using Whatshap".format(formatted_time))
        with concurrent.futures.ThreadPoolExecutor(THREADS) as executor:
            for chrome in CHR:
                command_WHATSHAP = [ WHATSHAP, "phase",
                            "--output", "{}/phased_{}.vcf.gz".format(PHASE_VCF_PATH,chrome),
                            "--reference", REFERENCE_FILE_PATH,
                            "--chromosome",chrome,
                            "--distrust-genotypes",
                            "--ignore-read-groups",
                            "{}/{}.vcf".format(PHASE_VCF_PATH,chrome),
                            BAM_FILE_PATH ]
                executor.submit(run_subprocess,command_WHATSHAP)
    with concurrent.futures.ThreadPoolExecutor(THREADS) as executor:
        for chrome in CHR:
            command_TABIX = ["tabix", "-f", "-p", "vcf", "{}/phased_{}.vcf.gz".format(PHASE_VCF_PATH,chrome)]
            executor.submit(run_subprocess,command_TABIX)

# Full alignment calling
#-----------------------------------------------------------------------------------------------------------------------
now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
print("[{}] STEP5 select candidates for full-alignment calling".format(formatted_time))
command_SelectQual = "gzip -fdc {}/pileup.vcf.gz | {} {} SelectQual --output_fn {} --var_pct_full {} --ref_pct_full {} --platform {} --vcf_fn {}".format(OUTPUT_FOLDER,PYPY,CLAIR3,CANDIDATE_BED_PATH,PRO,REF_PRO,PLATFORM,VCF_FILE_PATH)
subprocess.run(command_SelectQual,shell = True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)

with concurrent.futures.ThreadPoolExecutor(THREADS) as executor:
    for chrome in CHR:
        command_SelectCandidates = [PYPY, CLAIR3, "SelectCandidates",
                                    "--pileup_vcf_fn", "{}/pileup.vcf.gz".format(OUTPUT_FOLDER),
                                    "--split_folder", CANDIDATE_BED_PATH,
                                    "--ref_fn", REFERENCE_FILE_PATH,
                                    "--var_pct_full", PRO,
                                    "--ref_pct_full", REF_PRO,
                                    "--platform", PLATFORM,
                                    "--ctgName", chrome]
        executor.submit(run_subprocess,command_SelectCandidates)

now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
print("[{}] STEP6 call low-quality variants using full-alignment model".format(formatted_time))
file_num = subprocess.run("ls {}/FULL_ALN_FILE_* | wc -l".format(CANDIDATE_BED_PATH),capture_output=True, text=True,shell=True)
if int(file_num.stdout.strip()) == 0:
    print("[ERROR] No Candidate found! Exit in selecting full-alignment candidates")
    sys.exit()
else:
    subprocess.run("cat {}/FULL_ALN_FILE_* > {}/FULL_ALN_FILES".format(CANDIDATE_BED_PATH,CANDIDATE_BED_PATH), capture_output=True, text=True, shell=True)

with concurrent.futures.ThreadPoolExecutor(THREADS) as executor, open("{}/FULL_ALN_FILES".format(CANDIDATE_BED_PATH),"r") as FAFs:
    for line in FAFs:
        FAF = line.strip()
        command_CallVariantsFromCffi = [PYTHON, CLAIR3, "CallVariantsFromCffi",
                                        "--chkpnt_fn", FULL_ALIGNMENT_CHECKPOINT_PATH,
                                        "--bam_fn", BAM_FILE_PATH,
                                        "--call_fn", "{}/full_alignment_{}.vcf".format(FULL_ALIGNMENT_OUTPUT_PATH,FAF.split("/")[-1]),
                                        "--sampleName", SAMPLE,
                                        "--vcf_fn", VCF_FILE_PATH,
                                        "--ref_fn", REFERENCE_FILE_PATH,
                                        "--full_aln_regions", FAF,
                                        "--ctgName", FAF.split("/")[-1].split(".")[0],
                                        "--add_indel_length",
                                        "--no_phasing_for_fa", str(non_phase),
                                        "--minMQ", MIN_MQ,
                                        "--minCoverage", MIN_COV,
                                        "--phased_vcf_fn", "{}/phased_{}.vcf.gz".format(PHASE_VCF_PATH,FAF.split("/")[-1].split(".")[0]),
                                        "--gvcf", GVCF,
                                        "--enable_long_indel", ENABLE_LONG_INDEL,
                                        "--samtools", SAMTOOLS,
                                        "--use_gpu", USE_GPU,
                                        "--keep_iupac_bases", KEEP_IUPAC_BASES,
                                        "--platform", PLATFORM]
        executor.submit(run_subprocess,command_CallVariantsFromCffi)

command_Sortvcf = [PYPY, CLAIR3, "SortVcf",
    "--input_dir", FULL_ALIGNMENT_OUTPUT_PATH,
    "--vcf_fn_prefix", "full_alignment",
    "--output_fn", "{}/full_alignment.vcf".format(OUTPUT_FOLDER),
    "--sampleName", SAMPLE,
    "--ref_fn", REFERENCE_FILE_PATH,
    "--contigs_fn", "{}/CONTIGS".format(TMP_FILE_PATH),
    "--cmd_fn", CMD_FN]
subprocess.run(command_Sortvcf, capture_output=True, text=True)

full_alignment_vcf_num = subprocess.run("zcat {}/full_alignment.vcf.gz | grep -v '#' | wc -l".format(OUTPUT_FOLDER), capture_output=True, text=True, shell=True)
if int(full_alignment_vcf_num.stdout.strip()) == 0:
    print("[INFO] Exit in full-alignment variant calling")
    sys.exit()

if GVCF == "True":
    command_SortVcf = [PYPY, CLAIR3, "SortVcf",
                        "--input_dir",GVCF_TMP_PATH,
                        "--vcf_fn_suffix", ".tmp.gvcf",
                        "--output_fn", "{}/non_var.gvcf".format(GVCF_TMP_PATH),
                        "--ref_fn", REFERENCE_FILE_PATH,
                        "--contigs_fn", os.path.join(TMP_FILE_PATH,"CONTIGS"),
                        "--cmd_fn", CMD_FN]
    subprocess.run(command_Sortvcf, capture_output=True, text=True)

##Merge pileup and full alignment vcf
##-----------------------------------------------------------------------------------------------------------------------
now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
print("[{}] STEP7 Merge pileup VCF and full-alignment VCF".format(formatted_time))
SHOW_REF = str(args.print_ref_calls)
HAP_PRE = str(args.haploid_precise)
HAP_SEN = str(args.haploid_sensitive)
with concurrent.futures.ThreadPoolExecutor(THREADS) as executor:
    for chrome in CHR:
        command_MergeVcf = [PYPY, CLAIR3, "MergeVcf",
                            "--pileup_vcf_fn", "{}/pileup.vcf.gz".format(OUTPUT_FOLDER),
                            "--bed_fn_prefix", CANDIDATE_BED_PATH,
                            "--full_alignment_vcf_fn", "{}/full_alignment.vcf.gz".format(OUTPUT_FOLDER),
                            "--output_fn", "{}/merge_output/merge_{}.vcf".format(TMP_FILE_PATH,chrome),
                            "--platform", PLATFORM,
                            "--print_ref_calls", SHOW_REF,
                            "--gvcf", GVCF,
                            "--haploid_precise", HAP_PRE,
                            "--haploid_sensitive", HAP_SEN,
                            "--gvcf_fn", "{}/merge_output/merge_{}.gvcf".format(TMP_FILE_PATH,chrome),
                            "--non_var_gvcf_fn", "{}/non_var.gvcf".format(GVCF_TMP_PATH),
                            "--ref_fn", REFERENCE_FILE_PATH,
                            "--ctgName", chrome]
        executor.submit(run_subprocess,command_MergeVcf)

command_Sortvcf = [PYPY, CLAIR3, "SortVcf",
                    "--input_dir", "{}/merge_output".format(TMP_FILE_PATH),
                    "--vcf_fn_prefix", "merge",
                    "--output_fn", "{}/merge_output.vcf".format(OUTPUT_FOLDER),
                    "--sampleName", SAMPLE,
                    "--ref_fn", REFERENCE_FILE_PATH,
                    "--contigs_fn", "{}/CONTIGS".format(TMP_FILE_PATH),
                    "--cmd_fn", CMD_FN]
subprocess.run(command_Sortvcf,capture_output=True, text=True)

##Merge the output of genotype filling and DL
##-----------------------------------------------------------------------------------------------------------------------
now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
print("[{}] STEP8 Merge final result".format(formatted_time))
genotype_filling_vcf_list = {}
for file in os.listdir(Genotype_filling_fn):
    if len(file.split(".")) > 1:
        if file.split(".")[1] == "vcf" and file.split(".")[-1] == "gz":
            key = file.split(".")[0]
            if "chr" in file:
                if key == "chrX":
                    key = 23
                elif key == "chrY":
                    key = 24
                else:
                    key = int(key.replace("chr",""))
            else:
                if key == "X":
                    key = 23
                elif key == "Y":
                    key = 24
                else:
                    key = int(key)
            file = os.path.join(Genotype_filling_fn,file)
            genotype_filling_vcf_list[key] = file
genotype_filling_vcf_list = dict(sorted(genotype_filling_vcf_list.items()))
new_header_filter = pysam.VariantHeader()
new_header_info = pysam.VariantHeader()
new_header_format = pysam.VariantHeader()
new_header_contig = pysam.VariantHeader()
new_header = pysam.VariantHeader()
new_header.add_line('##fileformat=VCFv4.2')
new_header.add_line("##source=Genotype filled")
new_header.add_sample(args.sample_name)

for file in genotype_filling_vcf_list.values():
    vcftmp = pysam.VariantFile(file, "r")
    for record in vcftmp.header.records:
        if record.key == 'FILTER':
            if record not in new_header_filter.records:
                new_header_filter.add_record(record)
        elif record.key == 'INFO':
            if record not in new_header_info.records:
                new_header_info.add_record(record)
        elif record.key == 'FORMAT':
            if record not in new_header_format.records:
                new_header_format.add_record(record)
        elif record.key == 'contig':
            if record not in new_header_contig.records:
                new_header_contig.add_record(record)
    vcftmp.close()

for record in new_header_filter.records:
    new_header.add_record(record)
for record in new_header_info.records:
    new_header.add_record(record)
for record in new_header_format.records:
    new_header.add_record(record)
for record in new_header_contig.records:
    new_header.add_record(record)

Genotype_filling_final_vcf = os.path.join(OUTPUT_FOLDER,"Genotype_filling_final.vcf")
vcf_out = pysam.VariantFile(Genotype_filling_final_vcf, "w", header=new_header)
vcf_out.close()
with open(Genotype_filling_final_vcf,"a") as vcf_out:
    for i in genotype_filling_vcf_list.values():
        vcftmp = pysam.VariantFile(i, "r")
        for rec in vcftmp:
            vcf_out.write(str(rec))
        vcftmp.close()
subprocess.run("{} view {} -Oz -o {} ; {} -p vcf {} ; rm {}".format(BCFTOOLS,Genotype_filling_final_vcf,Genotype_filling_final_vcf+".gz",TABIX,Genotype_filling_final_vcf+".gz",Genotype_filling_final_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

final_merge = os.path.join(script_path,"Cyclone","Merge_vcf.py")
final_merge_vcf = os.path.join(OUTPUT_FOLDER,"Final_merged.vcf")
command_Merge = [PYTHON, final_merge,
                 "--vcf_qual_fn", "{}/merge_output.vcf.gz".format(OUTPUT_FOLDER),
                 "--vcf_gp_fn", Genotype_filling_final_vcf+".gz",
                 "--output", final_merge_vcf,
                 "--bcftools", BCFTOOLS,
                 "--tabix", TABIX]
subprocess.run(command_Merge, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
subprocess.run("{} sort {} -Oz -o {} -T {} ; {} -p vcf {} ; rm {}".format(BCFTOOLS,final_merge_vcf,final_merge_vcf+".gz",OUTPUT_FOLDER,TABIX,final_merge_vcf+".gz",final_merge_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
