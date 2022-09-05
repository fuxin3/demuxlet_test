#!/bin/bash

# ------------------------------------------------
# Shell script to generate dataset for demuxlet testing
# ------------------------------------------------

# to merge and index parsed BAM files

mkdir bam_merged

samtools merge bam_merged/bam_merged.bam bw_1.bam pbmc_1.bam pbmc_2.bam

samtools index bam_merged/bam_merged.bam

# to parse and merge cell barcode files

mkdir barcodes_merged

ln -s /SGRNJ06/randd/PROJECT/RD20102301_DZH/P22061404_YH_CART/20220801/C017T1-010-19-2FS/05.count/C017T1-010-19-2FS_filtered_feature_bc_matrix/barcodes.tsv  barcodes_merged/barcodes_bw_1.tsv
ln -s /SGRNJ06/randd/PROJECT/RD20102301_DZH/P22061404_YH_CART/20220801/C017T1-010-20FS/05.count/C017T1-010-20FS_filtered_feature_bc_matrix/barcodes.tsv barcodes_merged/barcodes_pbmc_1.tsv
ln -s /SGRNJ06/randd/PROJECT/RD20102301_DZH/P22061404_YH_CART/20220801/C017T1-010-21-2FS/05.count/C017T1-010-21-2FS_filtered_feature_bc_matrix/barcodes.tsv barcodes_merged/barcodes_pbmc_2.tsv

sed -i "s|\([A-Z]\+\)\-1|\1\-bw_1|g" barcodes_merged/barcodes_bw_1.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-pbmc_1|g" barcodes_merged/barcodes_pbmc_1.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-pbmc_2|g" barcodes_merged/barcodes_pbmc_2.tsv

cat barcodes_merged/barcodes_bw_1.tsv barcodes_merged/barcodes_pbmc_1.tsv barcodes_merged/barcodes_pbmc_2.tsv > barcodes_merged/barcodes_merged.tsv

#有重复
Rscript /generate_awk_lookup_tables_doublets.R 

#run doublets simulation
samtools view -h HGSOC/bam_merged/bam_merged.bam | \
awk \
'NR==1 { next } FNR==NR { a[$1]=$2; next } (i=gensub(/.*CB\:Z\:([A-Za-z]+\-[A-Za-z0-9]+).*/, "\\1", 1, $0)) in a { gsub(i, a[i]) }1' \
HGSOC/doublets_sims/20pc/lookup_table_doublets_HGSOC_20pc.tsv - | \
samtools view -bo HGSOC/doublets_sims/20pc/bam_merged_doublets_HGSOC_20pc.bam

#genotype
cellsnp-lite -s ../bam/bw_1.bam -b ../HGSOC/barcodes_merged/barcodes_bw_1.tsv -O ./bw_1 -p 10 --minMAF=0.01 --minCOUNT=50
cellsnp-lite -s ../bam/pbmc_1.bam -b ../HGSOC/barcodes_merged/barcodes_pbmc_1.tsv -O ./pbmc_1 -p 10 --minMAF=0.01 --minCOUNT=50
cellsnp-lite -s ../bam/pbmc_2.bam -b ../HGSOC/barcodes_merged/barcodes_pbmc_2.tsv -O ./pbmc_2 -p 10 --minMAF=0.01 --minCOUNT=50

gunzip -c ./bw_1/cellSNP.base.vcf.gz > ./bw_1/cellSNP.base-bgz.vcf
bgzip ./bw_1/cellSNP.base-bgz.vcf
gunzip -c ./pbmc_1/cellSNP.base.vcf.gz > ./pbmc_1/cellSNP.base-bgz.vcf
bgzip ./pbmc_1/cellSNP.base-bgz.vcf
gunzip -c ./pbmc_2/cellSNP.base.vcf.gz > ./pbmc_2/cellSNP.base-bgz.vcf
bgzip ./pbmc_2/cellSNP.base-bgz.vcf


# concatenate VCF files using vcftools (vcf-concat)

mkdir -p cellSNP_singlecell_merged

vcf-concat ./bw_1/cellSNP.base-bgz.vcf.gz ./pbmc_1/cellSNP.base-bgz.vcf.gz ./pbmc_2/cellSNP.base-bgz.vcf.gz > \
./cellSNP_singlecell_merged/cellSNP.base-merged.vcf


# run demuxlet
test_directory=/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/YiHe
/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/demuxlet/demuxlet \
--sam ${test_directory}/HGSOC/doublets_sims/20pc/bam_merged_doublets_HGSOC_20pc.bam
--vcf ${test_directory}/genotype/cellSNP_singlecell_merged/cellSNP.base-merged.vcf
--out ${test_directory}/output
--alpha 0.3
#--group-list 