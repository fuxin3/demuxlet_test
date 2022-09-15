# demuxlet

## Introduction

A tool for Multiplexed droplet detection in single-cell RNA-sequencing using natural genetic variation.

## Env

`conda activate demuxlet`

## Data prepare

- bam file, must be sorted by coordinates and indexed.
- vcf file, containing the individual genotypes (GT), posterior probability (GP), or genotype likelihood (PL).

## Usage

```
    cd /SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/data

    /SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/demuxlet/demuxlet \
    --sam doublets_sims/20pc/bam_merged_doublets_20pc.bam \
    --vcf genotype/cellSNP.base-merged.vcf \
    --group-list barcodes_merged/barcodes_merged.tsv \
    --write-pair \
    --out output3 \
    --alpha 0.5 
```

## Test

directory: /SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/data

script: [work.sh](https://github.com/fuxin3/demuxlet_test/blob/main/work.sh) 

## Parameters&Output detail for demuxlet

refer to [demuxlet](https://github.com/statgen/demuxlet)

