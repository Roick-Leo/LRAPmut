# LRAPmut
LRAPMut is a germline small variant caller for long-reads. 

# Install
```
conda create -n LRAPmut
conda install -c bioconda bcftools samtools clair3 python==3.9.0 -y
git clone https://github.com/Roick-Leo/LRAPmut.git
chmod -R a+x ./LRAPmut
```

# Usage
## General Usage
```
./lrapmut.bin \
    --env_dir=${The_bin_dir_of_LRAPmut_conda_env}
    --bam_fn=${BAM} \
    --ref_fn=${REF} \
    --threads=${THREADS} \  		     
    --platform="ont" \  
    --model_path=${MODEL_PREFIX} \
    --output=${OUTPUT_DIR} \
    --reference_panel_fn=${REF_PANEL_DIR} \
    --reference_map_fn=${REF_MAP_DIR} \
```

## Demo
```
minimap2 -ax map-ont --secondary=no -t 15 /path/to/ref /path/to/fastq | samtools view -@ 15 -bS | samtools sort -@ 15 -o /path/to/output_sorted.bam

/path/to/LRAPmut.bin --env_dir /path/to/LRAPmut_conda_env_dir/bin --bam_fn /path/to/bam_file --ref_fn /path/to/GRch38_ref.fa --model_path /path/to/model_dir --platform "ont" --threads 15 --output /path/to/output_dir --reference_panel_fn /path/to/refrence_panel_dir --reference_map_fn /path/to/genetic_maps.b38
```

## Options
### Required parameters:
```
  --env_dir                 The bin dir of LRAPmut conda env
  --bam_fn=FILE             BAM file input. The input file must be samtools indexed.
  --ref_fn=FILE             FASTA reference file input. The input file must be samtools indexed.
  --model_path=STR          The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00002, pileup.data-00001-of-00002 pileup.index, full_alignment.data-00000-of-00002, full_alignment.data-00001-of-00002  and full_alignment.index).
  --threads=INT             Max threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The chunks being processed simultaneously is ceil($threads/4)*3. 3 is the overloading factor.
  --platform=STR            Select the sequencing platform of the input.
  --output=PATH             VCF/GVCF output directory.
```

# How to get ref-panel
## Download reference panel
```
mkdir -p /path/to/download_dir
cd /path/to/download_dir
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr{1..22}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz{,.tbi}
```
## Remove NA12878 (and family) from the reference pane
```
conda activate LRAPmut
mkdir -p /path/to/reference_panel_dir
./home/helei/biosoft/LRAPmut/Greate_ref_panel.bin \
  --input_dir /path/to/ref_panel_download_dir \
  --output_dir /path/to/reference_panel_dir \
  --sample NA12878,NA12891,NA12892 \
```

# How to get genetic maps
```
mkdir -p /path/to/genetic_maps_dir
cd /path/to/genetic_maps_dir
wget -c https://raw.githubusercontent.com/odelaneau/shapeit4/master/maps/genetic_maps.b38.tar.gz
tar -xzvf file.tar.gz
```