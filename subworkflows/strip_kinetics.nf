#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// inputs
params.sample_sheet = ''

// output directory
params.outdir = "${workflow.launchDir}/fiberseq_output"

process strip_kinetics {
    publishDir "$params.outdir/stripped_kinetics", mode: 'copy'
    cpus 16
    memory '64 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sample_id), path(input_bam), val(ref_name)
    output:
    tuple val(sample_id), path("${input_bam.baseName}.stripped.bam")
    script:
    """
    samtools view -b -@7 \
    --remove-tag fi,ri,fp,rp,ip,pw \
    ${input_bam} > ${input_bam.baseName}.stripped.bam
    """
}
workflow{    
 
    samplesheet_ch = channel.fromPath("${params.sample_sheet}")
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> tuple(row[0], file(row[1]), row[2])}

    // samplesheet_ch.view()        
    strip_kinetics(samplesheet_ch)
}