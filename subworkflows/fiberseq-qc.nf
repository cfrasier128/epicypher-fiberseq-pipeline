#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process create_qc_table {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}_qc_table.tbl"), emit: qc_table
    env ("tech"), emit: tech
    env ("nreads"), emit: nreads

    script:
    """
    tech=`samtools view -H ${bam_file} | grep -e "PL:" | awk 'NR == 1' | tr '\t' '\n' | awk '\$1 ~ /^PL:/' | cut -f2 -d':' | tr '[[:upper:]]' '[[:lower:]]' || true`
    nreads=`samtools view -c ${bam_file}`
    export tech
    export nreads
    ft extract -t 8 ${bam_file} --all - \
    | /opt/fiberseq/details/cutnm 5mC,ec,fiber,fiber_length,m6a,msp_lengths,msp_starts,nuc_lengths,nuc_starts,rq,total_5mC_bp,total_AT_bp,total_m6a_bp \
    > ${sample_id}_qc_table.tbl
    """
}

process plot_msp_lengths {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(qc_table)

    output:
    tuple val(sample_id), path("${sample_id}.msp_lengths.pdf"), path("${sample_id}.msp_lengths.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-msp-lengths.sh ${sample_id} ${qc_table} ${sample_id}.msp_lengths.pdf ${sample_id}.msp_lengths.intermediate.stat.txt
    """
}

process plot_nuc_lengths {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(qc_table)

    output:
    tuple val(sample_id), path("${sample_id}.nuc_lengths.pdf"), path("${sample_id}.nuc_lengths.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-nuc-lengths.sh ${sample_id} ${qc_table} ${sample_id}.nuc_lengths.pdf ${sample_id}.nuc_lengths.intermediate.stat.txt
    """
}

process plot_6ma {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(qc_table)

    output:
    tuple val(sample_id), path("${sample_id}.m6a_per_read.pdf"), path("${sample_id}.m6a_per_read.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-number-m6a-per-read.sh ${sample_id} ${qc_table} ${sample_id}.m6a_per_read.pdf ${sample_id}.m6a_per_read.intermediate.stat.txt ${sample_id}.ccs_passes.pdf ${sample_id}.ccs_passes.intermediate.stat.txt
    """
}

process plot_nucs_per_read {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(qc_table)

    output:
    tuple val(sample_id), path("${sample_id}.number_nucs_per_read.pdf"), path("${sample_id}.number_nucs_per_read.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-number-nucs-per-read.sh ${sample_id} ${qc_table} ${sample_id}.number_nucs_per_read.pdf ${sample_id}.number_nucs_per_read.intermediate.stat.txt
    """
}

process plot_5mCs_per_read {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(qc_table)

    output:
    tuple val(sample_id), path("${sample_id}.number_cpgs_per_read.pdf"), path("${sample_id}.number_cpgs_per_read.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-number-mcpgs-per-read.sh ${sample_id} ${qc_table} ${sample_id}.number_cpgs_per_read.pdf ${sample_id}.number_cpgs_per_read.intermediate.stat.txt
    """
}

process plot_readlength_per_nuc {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(qc_table)

    output:
    tuple val(sample_id), path("${sample_id}.readlength_per_nuc.pdf"), path("${sample_id}.readlength_per_nuc.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-readlength-per-nuc.sh ${sample_id} ${qc_table} ${sample_id}.readlength_per_nuc.pdf ${sample_id}.readlength_per_nuc.intermediate.stat.txt
    """
}

process plot_readlengths {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(qc_table)
    env "tech"

    output:
    tuple val(sample_id), path("${sample_id}.readlengths.pdf"), path("${sample_id}.readlengths.intermediate.stat.txt")

    script:
    """
    #if [[ "\${tech}" == "pacbio" ]]; then
    #    export max_scale=25000
    #fi
    #if [[ "\${tech}" == "ont" ]]; then
        export max_scale=50000
    #fi
    /opt/fiberseq/details/make-plot-readlengths.sh ${sample_id} ${qc_table} \${max_scale} ${sample_id}.readlengths.pdf ${sample_id}.readlengths.intermediate.stat.txt
    """
}

process plot_msp_resolution {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(qc_table)

    output:
    tuple val(sample_id), path("${sample_id}.msp_resolution.pdf"), path("${sample_id}.msp_resolution.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-msp-resolution.sh ${sample_id} ${qc_table} ${sample_id}.msp_resolution.pdf ${sample_id}.msp_resolution.intermediate.stat.txt
    """
}

process plot_read_quality {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(qc_table)

    output:
    tuple val(sample_id), path("${sample_id}.readquality.pdf"), path("${sample_id}.readquality.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-rq.sh ${sample_id} ${qc_table} ${sample_id}.readquality.pdf ${sample_id}.readquality.intermediate.stat.txt
    """
}

process plot_autocorrelation {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(bam_file)
    env "tech"

    output:
    tuple val(sample_id), path("${sample_id}.autocorrelation.pdf"), path("${sample_id}.autocorrelation.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-autocorrelation.v2.sh ${sample_id} ${bam_file} ${sample_id}.autocorrelation.pdf ${sample_id}.autocorrelation.intermediate.stat.txt \$tech
    """
}

process plot_randfibers {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.randfibers.pdf"), path("${sample_id}.randfibers.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-rand-fibers.sh ${sample_id} ${bam_file} 0 20000 ${sample_id}.randfibers.pdf ${sample_id}.randfibers.intermediate.stat.txt
    """
}

process plot_zoomed_randfibers {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.randfibers.2K-4K.pdf"), path("${sample_id}.randfibers.2K-4K.intermediate.stat.txt")

    script:
    """
    /opt/fiberseq/details/make-plot-rand-fibers.sh ${sample_id} ${bam_file} 2000 4000 ${sample_id}.randfibers.2K-4K.pdf ${sample_id}.randfibers.2K-4K.intermediate.stat.txt
    """
}

process join_qc {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(msp_pdf), path(msp_txt), path(nuc_pdf), path(nuc_txt), path(ma_pdf), path(ma_txt), path(nucs_per_read_pdf), path(nucs_per_read_txt), path(cpgs_per_read_pdf), path(cpgs_per_read_txt), path(readlength_per_nuc_pdf), path(readlength_per_nuc_txt), path(readlengths_pdf), path(readlengths_txt), path(msp_resolution_pdf), path(msp_resolution_txt), path(autocorrelation_pdf), path(autocorrelation_txt), path(randfibers_pdf), path(randfibers_txt), path(zoomed_randfibers_pdf), path(zoomed_randfibers_txt)
    env "nreads"

    output:
    path "${sample_id}.qc_stats.txt"

    script:
    """
    echo \$nreads \
  | awk '{ printf "# Note: ***Unaligned Reads***\n# Stats: nreads\nNumber(Reads)=%s\n\n", \$1; }' \
  | cat - *txt \
 > ${sample_id}.qc_stats.txt
    """
}

process make_html {
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fiberseq-qc:latest'

    input:
    tuple val(sample_id), path(msp_pdf), path(msp_txt), path(nuc_pdf), path(nuc_txt), path(ma_pdf), path(ma_txt), path(nucs_per_read_pdf), path(nucs_per_read_txt), path(cpgs_per_read_pdf), path(cpgs_per_read_txt), path(readlength_per_nuc_pdf), path(readlength_per_nuc_txt), path(readlengths_pdf), path(readlengths_txt), path(msp_resolution_pdf), path(msp_resolution_txt), path(autocorrelation_pdf), path(autocorrelation_txt), path(randfibers_pdf), path(randfibers_txt), path(zoomed_randfibers_pdf), path(zoomed_randfibers_txt)
    env "nreads"

    output:
    path "${sample_id}.qc.html"
    path "${sample_id}.overview.html"
    path "${sample_id}*.png"

    script:
    """
    /opt/fiberseq/details/make-html.tcsh \
    ${sample_id} \
    \$nreads \
    ${sample_id}.overview.html \
    ${sample_id}.qc.html \
    
    """
}

workflow fiberseq_qc_workflow {
    take:
    msp_bams

    main:
    input_ch = msp_bams.map { row -> tuple(row[0], row[1]) }
    create_qc_table(input_ch)
    plot_msp_lengths(create_qc_table.out.qc_table)
    plot_nuc_lengths(create_qc_table.out.qc_table)
    plot_6ma(create_qc_table.out.qc_table)
    plot_nucs_per_read(create_qc_table.out.qc_table)
    plot_5mCs_per_read(create_qc_table.out.qc_table)
    plot_readlength_per_nuc(create_qc_table.out.qc_table)
    plot_readlengths(create_qc_table.out.qc_table, create_qc_table.out.tech)
    plot_msp_resolution(create_qc_table.out.qc_table)
    //plot_read_quality(create_qc_table.out.qc_table)
    plot_autocorrelation(input_ch, create_qc_table.out.tech)
    plot_randfibers(input_ch)
    plot_zoomed_randfibers(input_ch)

    all_qc_file_ch = plot_msp_lengths.out
        .combine(
            plot_nuc_lengths.out,
            by: 0
        )
        .combine(
            plot_6ma.out,
            by: 0
        )
        .combine(
            plot_nucs_per_read.out,
            by: 0
        )
        .combine(
            plot_5mCs_per_read.out,
            by: 0
        )
        .combine(
            plot_readlength_per_nuc.out,
            by: 0
        )
        .combine(
            plot_readlengths.out,
            by: 0
        )
        .combine(
            plot_msp_resolution.out,
            by: 0
        )
        .combine(
            plot_autocorrelation.out,
            by: 0
        )
        .combine(
            plot_randfibers.out,
            by: 0
        )
        .combine(
            plot_zoomed_randfibers.out,
            by: 0
        )

    join_qc(all_qc_file_ch, create_qc_table.out.nreads)
    make_html(all_qc_file_ch, create_qc_table.out.nreads)
}
