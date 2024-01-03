# LRAPmut
LRAPMut is a germline small variant caller for long-reads. 

# Usage
## General Usage
```
python ./lrapmut.py \
    --bam_fn=${BAM} \
    --ref_fn=${REF} \
    --threads=${THREADS} \  		     
    --platform="ont" \  
    --model_path=${MODEL_PREFIX} \
    --output=${OUTPUT_DIR} \
    --reference_panel_fn=${REF_PANEL_DIR} \
    --reference_map_fn=${REF_MAP_DIR} \
    --tools=${TOOLS_DIR}
```

## Options
### Required parameters:
```
  -b, --bam_fn=FILE             BAM file input. The input file must be samtools indexed.
  -f, --ref_fn=FILE             FASTA reference file input. The input file must be samtools indexed.
  -m, --model_path=STR          The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00002, pileup.data-00001-of-00002 pileup.index, full_alignment.data-00000-of-00002, full_alignment.data-00001-of-00002  and full_alignment.index).
  -t, --threads=INT             Max threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The chunks being processed simultaneously is ceil($threads/4)*3. 3 is the overloading factor.
  -p, --platform=STR            Select the sequencing platform of the input.
  -o, --output=PATH             VCF/GVCF output directory.
```