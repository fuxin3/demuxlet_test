#!/bin/bash

# ------------------------------------------------
#script to generate dataset for demuxlet testing
# ------------------------------------------------


test_directory=/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/YiHe_sub
conda activate demuxlet

#取chr21 的数据来测试，减少数据运行时间
cd ${test_directory}
mkdir bam

samtools view  -hb ../YiHe/bam/bw_1.bam 21 > bam/bw_1.chr21.bam
samtools view  -hb ../YiHe/bam/pbmc_1.bam 21 > bam/pbmc_1.chr21.bam
samtools view  -hb ../YiHe/bam/pbmc_2.bam 21 > bam/pbmc_2.chr21.bam

samtools index bam/bw_1.chr21.bam
samtools index bam/pbmc_1.chr21.bam
samtools index bam/pbmc_2.chr21.bam

# to merge and index parsed BAM files

mkdir bam_merged

samtools merge bam_merged/bam_merged.bam bam/bw_1.chr21.bam bam/pbmc_1.chr21.bam bam/pbmc_2.chr21.bam

samtools index bam_merged/bam_merged.bam

# to parse and merge cell barcode files
mkdir barcodes_merged

## run celescope
celescope rna count --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --outdir celescope/bw_1/05.count --sample bw_1 --bam bam/bw_1.chr21.bam
celescope rna count --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --outdir celescope/pbmc_1/05.count --sample pbmc_1 --bam bam/pbmc_1.chr21.bam
celescope rna count --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --outdir celescope/pbmc_2/05.count --sample pbmc_2 --bam bam/pbmc_2.chr21.bam


ln -s /SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/YiHe_sub/celescope/bw_1/05.count/bw_1_matrix_10X/barcodes.tsv  barcodes_merged/barcodes_bw_1.tsv
ln -s /SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/YiHe_sub/celescope/pbmc_1/05.count/pbmc_1_matrix_10X/barcodes.tsv barcodes_merged/barcodes_pbmc_1.tsv
ln -s /SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/YiHe_sub/celescope/pbmc_2/05.count/pbmc_2_matrix_10X/barcodes.tsv barcodes_merged/barcodes_pbmc_2.tsv

sed -i "s|\([A-Z]\+\)\-1|\1\-bw_1|g" barcodes_merged/barcodes_bw_1.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-pbmc_1|g" barcodes_merged/barcodes_pbmc_1.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-pbmc_2|g" barcodes_merged/barcodes_pbmc_2.tsv

cat barcodes_merged/barcodes_bw_1.tsv barcodes_merged/barcodes_pbmc_1.tsv barcodes_merged/barcodes_pbmc_2.tsv > barcodes_merged/barcodes_merged.tsv

# 
Rscript /generate_awk_lookup_tables_doublets.R 
# 

##run doublets simulation
samtools view -h bam_merged/bam_merged.bam | \
awk \
'NR==1 { next } FNR==NR { a[$1]=$2; next } (i=gensub(/.*CB\:Z\:([A-Za-z]+\-[A-Za-z0-9]+).*/, "\\1", 1, $0)) in a { gsub(i, a[i]) }1' \
doublets_sims/20pc/lookup_table_doublets_chr21_20pc.tsv - | \
samtools view -bo doublets_sims/20pc/bam_merged_doublets_chr21_20pc.bam

#genotype
#cellsnp-lite call不出变异 
#改用bcftools
samtools mpileup -t DP -t SP -t DP4 -uf \
/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa \
bam/bw_1.chr21.bam | \
bcftools call -f GQ -f GP -mv > genotype/bw_1.vcf

samtools mpileup -t DP -t SP -t DP4 -uf \
/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa bam/pbmc_1.chr21.bam | \
bcftools call -f GQ -f GP -mv > genotype/pbmc_1.vcf

samtools mpileup -t DP -t SP -t DP4 -uf \
/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa bam/pbmc_2.chr21.bam | \
bcftools call -f GQ -f GP -mv > genotype/pbmc_2.vcf


#vcf文件的合并
#vcf-concat ./genotype/bw_1.1.vcf ./genotype/pbmc_1.1.vcf ./genotype/pbmc_2.1.vcf > ./genotype/cellSNP.base-merged.vcf
bgzip bw_1.1.vcf
bgzip pbmc_1.1.vcf
bgzip pbmc_2.1.vcf

bcftools index -t bw_1.1.vcf.gz
bcftools index -t pbmc_1.1.vcf.gz
bcftools index -t pbmc_2.1.vcf.gz

vcf-merge ./genotype/bw_1.1.vcf.gz ./genotype/pbmc_1.1.vcf.gz ./genotype/pbmc_2.1.vcf.gz > ./genotype/cellSNP.base-merged.vcf



#run demuxlet
/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/demuxlet/demuxlet \
--sam doublets_sims/20pc/bam_merged_doublets_chr21_20pc.bam
--vcf genotype/cellSNP.base-merged.vcf
--out output
--alpha 0.3