#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process create_pileups {
    label 'large'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index)

    output:
    tuple val(samp_name), path("*.tsv.gz"), val(ref_name), emit: pileups

    script:
    """
    ft pileup \
        --m6a \
        --cpg \
        -t ${task.cpus} \
        --ftx "len(msp)>10" \
        ${aligned_bam} \
    | awk -v OFS="\t" -v FS="\t" '{print \$1,\$2,\$3,\$9/(\$4+0.1),\$10/(\$4+0.1),\$7/(\$4+0.1)}' \
    | gzip -c > ${samp_name}.pileup_all.tsv.gz
    """
}

process create_5mC_pileup_methbat {
    label 'large'
    container 'quay.io/pacbio/methbat:1.1.0_build2'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index), path(ref_fasta), path(ref_fai)

    output:
    tuple val(samp_name), path("*.tsv.gz"), val(ref_name), emit: pileups

    script:
    """
    aligned_bam_to_cpg_scores \
      --threads ${task.cpus} \
      --bam ${aligned_bam} \
      --output-prefix ${samp_name}.${ref_name} \
      --min-mapq 1 \
      --min-coverage 4 \
      --skip-6mA
      --skip-5hmC
      --pileup-mode count

    gunzip -c ${samp_name}.${ref_name}.5mC.bed.gz \
    | awk -v OFS=\$'\t' '!/^#/ {print \$1, \$2, \$3, \$7/100}' \
    | gzip -c > ${samp_name}.cpgpileup.tsv.gz
    """
}

process pileupbedgraphtobigwig {
    publishDir "${params.outdir}/4_pileups/${samp_name}", mode: 'copy'
    label 'medium'
    container 'quay.io/pacbio/bigtools:3844b58_build1'

    input:
    tuple val(samp_name), path(bedgraph), val(ref_name), path(ref_fai), val(feature), val(col_num)

    output:
    tuple val(samp_name), path("*.bw")

    script:
    """
    cut -f 1,2 ${ref_fai} > chromsizes
    zcat ${bedgraph} \
    | cut -f 1,2,3,${col_num} \
    | grep -v '^#' \
    | sort -k1,1 -k2,2n \
    > temp.bedgraph
    
    bedgraphtobigwig \
        --nthreads ${task.cpus} \
        temp.bedgraph \
        chromsizes \
        ${samp_name}.${feature}.bw
    """
}

workflow create_bigwigs {
    take:
    msp_bams
    references_ch

    main:
        // If --create_bigwigs is set in command line, create pileups and bigwigs
        create_pileups(msp_bams)
        // -> samp_name, pileups, ref_name
        // create a bedgraph of 6ma calling and calculate percent 6ma coverage for each base
        create_5mC_pileup_methbat(msp_bams.combine(references_ch, by: 2).map { row -> tuple(row[1], row[2], row[0], row[3], row[4], row[5]) })
        create_pileups.out.pileups
            .combine(references_ch, by: 2)
            .map { row -> tuple(row[1], row[2], row[0], row[4]) }
            .set { pileup_withref_ch }
        // -> samp_name, pileups, ref_name, ref_fai
        // convert the pileup bedgraph into a bigwig for downstream purposes
        create_5mC_pileup_methbat.out.pileups
            .combine(references_ch, by: 2)
            .map { row -> tuple(row[1], row[2], row[0], row[4], "perccpg", 4) }
            .set { cpg_tools_pileup_bedgraph_ch }
        pileup_withref_ch
            .map { row -> tuple(row[0], row[1], row[2], row[3], "perc6ma", 4) }
            .concat(pileup_withref_ch.map { row -> tuple(row[0], row[1], row[2], row[3], "percnuc", 6) })
            .concat(cpg_tools_pileup_bedgraph_ch)
            .set { pileup_bedgraph_ch }
        pileupbedgraphtobigwig(pileup_bedgraph_ch)
}