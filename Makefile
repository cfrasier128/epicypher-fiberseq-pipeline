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

input_template:
	mkdir -p inputs
	cat > inputs/sample_sheet_template.tsv << EOL
samp_name\tbam_path\tref_name
Sample1\t/data/Sample1.run1.bam\thg38
Sample1\t/data/Sample1.run2.bam\thg38
Sample2\t/data/Sample2.bam\thg38
EOL

inputs: references input_template

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