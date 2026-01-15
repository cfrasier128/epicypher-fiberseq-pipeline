#!/usr/bin/env bash
set -euo pipefail

# Human GRCh38_no_alt_analysis_set with decoys and specific gene masks
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt
mkdir -p references/hg38
cd references/hg38 || exit 1
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
gunzip -c GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz > hg38.fasta
samtools faidx hg38.fasta
cd ../../

# Human T2T-CHM13 v2.0 with rCRS mitochondrial sequence and masked chrY PARs
# https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/README.txt
mkdir -p references/chm13
cd references/chm13 || exit 1
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz
gunzip -c chm13v2.0_maskedY_rCRS.fa.gz > chm13.fasta
samtools faidx chm13.fasta
cd ../../

# Mouse GRCm39 assembly analysis set
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/seqs_for_alignment_pipelines/README_analysis_sets.txt
mkdir -p references/mm39
cd references/mm39 || exit 1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/seqs_for_alignment_pipelines/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz
gunzip -c GCA_000001635.9_GRCm39_full_analysis_set.fna.gz > mm39.fasta
samtools faidx mm39.fasta
cd ../../

# Mouse mm10 assembly initial release
# https://hgdownload.gi.ucsc.edu/goldenPath/mm10/bigZips/
mkdir -p references/mm10
cd references/mm10 || exit 1
wget https://hgdownload.gi.ucsc.edu/goldenPath/mm10/bigZips/initial/mm10.fa.gz
gunzip -c mm10.fa.gz > mm10.fasta
samtools faidx mm10.fasta
cd ../../
