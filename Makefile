SHELL:=/bin/bash
export NXF_VER:=25.10.3
export PATH:=$(CURDIR):$(PATH)

export SINGULARITY_CACHEDIR:=$(CURDIR)/singularity_cache
export NXF_SINGULARITY_CACHEDIR:=$(SINGULARITY_CACHEDIR)

install: ./nextflow

./nextflow:
	curl -fsSL get.nextflow.io | bash

update:
	nextflow self-update

references:
	bash ./prepare_references.sh

check:
	nextflow lint main.nf
	nextflow lint subworkflows/fiberseq-qc.nf

format:
	nextflow lint -format main.nf
	nextflow lint -format subworkflows/fiberseq-qc.nf

clean:
	rm -vf .nextflow.log*
	rm -vrf .nextflow*
	rm -vf dag-*
	rm -vf report-*
	rm -vf timeline-*
	rm -vf trace-*
	rm -vrf work/*
	rm -vrf results/*

clean-singularity-cache:
	rm -vf singularity_cache/*.img

# nextflow run main.nf --sample_sheet inputs/sample_sheet.tsv -profile debug,slurm --pb_qc --phase_reads --create_bigwigs