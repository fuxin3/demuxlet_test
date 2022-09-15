#!/bin/bash

test_directory=/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/     
conda activate demuxlet

cd test_directory

# download test data
curl ftp.ebi.ac.uk/vol1/fastq/ERR386/002/ERR3863702/ERR3863702_1.fastq.gz -o ./data/souporcell/fastq/babz3_001/babz3_S1_L001_R1_001.fastq.gz
curl ftp.ebi.ac.uk/vol1/fastq/ERR386/002/ERR3863702/ERR3863702_2.fastq.gz -o ./data/souporcell/fastq/babz3_001/babz3_S1_L001_R2_001.fastq.gz
curl ftp.ebi.ac.uk/vol1/fastq/ERR386/003/ERR3863703/ERR3863703_1.fastq.gz -o ./data/souporcell/fastq/babz3_002/babz3_S1_L001_R1_002.fastq.gz
curl ftp.ebi.ac.uk/vol1/fastq/ERR386/003/ERR3863703/ERR3863703_2.fastq.gz -o ./data/souporcell/fastq/babz3_002/babz3_S1_L001_R2_002.fastq.gz

# run cellranger
cellranger count --id=babz3_001 \
--transcriptome=/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/data/Homo_sapiens.GRCh38.92 \
--fastqs=/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/data/souporcell/fastq/babz3_001 \
----sample=babz3 \
--nosecondary \
--jobmode=local \
--localcores=10 \
--localmem=50 \
--localvmem=100

cellranger count --id=babz3_002 \
--transcriptome=/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/data/Homo_sapiens.GRCh38.92 \
--fastqs=/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/data/souporcell/fastq/babz3_002 \
----sample=babz3 \
--nosecondary \
--jobmode=local \
--localcores=10 \
--localmem=50 \
--localvmem=100


# fix bam barcode
python /SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/fix_BC_for_bamfile.py \
data/babz3_001/outs/possorted_genome_bam.bam \
data/bam/babz3_001.fixBC.bam \
babz3_001

python /SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/fix_BC_for_bamfile.py \
data/babz3_002/outs/possorted_genome_bam.bam \
data/bam/babz3_001.fixBC.bam \
babz3_002

# to merge and index parsed BAM files

mkdir bam_merged

samtools merge data/bam_merged/bam_merged.bam data/bam/babz3_001.fixBC.bam data/bam/babz3_002.fixBC.bam

samtools index data/bam_merged/bam_merged.bam

# to parse and merge cell barcode files
mkdir data/barcodes_merged

cp data/babz3_001/outs/filtered_feature_bc_matrix/barcodes.tsv.gz data/barcodes_merged/babz3_001_barcodes.tsv.gz
cp data/babz3_002/outs/filtered_feature_bc_matrix/barcodes.tsv.gz data/barcodes_merged/babz3_002_barcodes.tsv.gz

gunzip data/barcodes_merged/babz3_001_barcodes.tsv.gz
gunzip data/barcodes_merged/babz3_002_barcodes.tsv.gz

sed -i "s|\([A-Z]\+\)\-1|\1\-babz3_001|g" data/barcodes_merged/babz3_001_barcodes.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-babz3_002|g" data/barcodes_merged/babz3_001_barcodes.tsv

cat data/barcodes_merged/babz3_001_barcodes.tsv data/barcodes_merged/babz3_002_barcodes.tsv > data/barcodes_merged/barcodes_merged.tsv  

# run doublets simulation

Rscript /SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/test/generate_awk_lookup_tables_doublets.R

samtools view -h bam_merged/bam_merged.bam | \
awk \
'NR==1 { next } FNR==NR { a[$1]=$2; next } (i=gensub(/.*CB\:Z\:([A-Za-z]+\-[A-Za-z0-9]+).*/, "\\1", 1, $0)) in a { gsub(i, a[i]) }1' \
doublets_sims/20pc/lookup_table_doublets_20pc.tsv - | \
samtools view -bo doublets_sims/20pc/bam_merged_doublets_20pc.bam

# genotype calling with bcftools
samtools mpileup -t DP -t SP -t DP4 -uf \
/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa data/bam/babz3_001.fixBC.bam | \
bcftools call -f GQ -f GP -mv > data/genotype/babz3_001.vcf

samtools mpileup -t DP -t SP -t DP4 -uf \
/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa data/bam/babz3_002.fixBC.bam | \
bcftools call -f GQ -f GP -mv > data/genotype/babz3_002.vcf

bgzip data/genotype/babz3_001.vcf
bgzip data/genotype/babz3_002.vcf

bcftools index -t data/genotype/babz3_001.vcf.gz
bcftools index -t data/genotype/babz3_002.vcf.gz

vcf-merge data/genotype/data/genotype/babz3_001.vcf.gz data/genotype/babz3_002.vcf.gz > data/genotype/cellSNP.base-merged.vcf

#run demuxlet
/SGRNJ03/randd/user/fuxin/PROJECTS/demultiplex/demuxlet/demuxlet \
--sam doublets_sims/20pc/bam_merged_doublets_20pc.bam
--vcf genotype/cellSNP.base-merged.vcf
--group-list barcodes_merged/barcodes_merged.tsv
--write-pair
--out output
--alpha 0.3
