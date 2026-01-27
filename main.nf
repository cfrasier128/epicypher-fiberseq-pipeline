#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// inputs
params.sample_sheet = ''
params.ref_path = "${workflow.launchDir}/references"

// default parameters
params.ref_name = 't2t'
params.confidence_ml_val = '250'
params.minimum_msp_dist = '10'

// optional steps
params.pb_qc = false
params.phase_reads = false
params.create_bigwigs = false
params.debug = false

// output directory
params.outdir = "${workflow.launchDir}/results"

// Grab subworkflows
include { fiberseq_qc_workflow } from './subworkflows/fiberseq-qc.nf'

process align_bams {
    label 'very_large'
    container 'quay.io/pacbio/pbmm2:1.17.0_build1'

    input:
    tuple val(samp_name), path(input_bam), val(ref_name), path(ref_fasta)

    output:
    tuple val(samp_name), path("*aligned.bam"), val(ref_name), path("*aligned.bam.bai"), emit: aligned_bam

    script:
    """
    out_prefix=\$(basename ${input_bam} .bam)

    pbmm2 align \
        --preset HIFI --log-level INFO \
        --strip \
        --num-threads ${task.cpus} \
        --sort --sort-memory 4G --bam-index BAI \
        --sample ${samp_name} \
        ${ref_fasta} \
        ${input_bam} \
        \${out_prefix}.aligned.bam
    """
}

process merge_bams {
    label 'large'
    container 'cfrasier/epi-pacbio:latest'

    input:
    tuple val(samp_name), path(input_bam_files), val(ref_name), val(input_bam_index_files)

    output:
    tuple val(samp_name), path("${samp_name}.${ref_name}.aligned.bam"), val(ref_name), path("${samp_name}.${ref_name}.aligned.bam.bai"), emit: merged_bams

    script:
    """
    pbmerge \
        --no-pbi --num-threads ${task.cpus} \
        -o ${samp_name}.${ref_name}.aligned.bam \
        ${input_bam_files.join(' ')}
    samtools index -@ ${task.cpus - 1} ${samp_name}.${ref_name}.aligned.bam
    """
}

process pacbio_qc {
    publishDir "${params.outdir}/sequencing_qc/${samp_name}", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-pacbio:latest'

    input:
    tuple val(samp_name), path(raw_bam), val(ref_name), path(bam_index)

    output:
    tuple path("ccs.report.json"), path("*plot.png"), emit: pacbio_qc_reports

    script:
    """
    pbindex ${raw_bam} -j ${task.cpus};
    dataset create --type ConsensusReadSet --name \$(echo ${samp_name} | cut -f 2 -d ".") ${samp_name}.hifi_reads.ccsreadset.xml ${raw_bam};
    runqc-reports -b --pdf-report ${samp_name}.pbqc.report.pdf ${samp_name}.hifi_reads.ccsreadset.xml
    """
}

process deepvariant {
    label 'very_large'
    container 'google/deepvariant:1.9.0'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index), path(ref_fasta), path(ref_fai)

    output:
    tuple val(samp_name), path("*vcf.gz"), val(ref_name), path("*vcf.gz.tbi"), emit: vcfs

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${ref_fasta} \
        --reads ${aligned_bam} \
        --output_vcf ${samp_name}.${ref_name}.deepvariant.vcf.gz \
        --num_shards ${task.cpus}
    """
}

process pbsv {
    label 'large'
    container 'quay.io/pacbio/pbsv:2.11.0_build1'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index), path(ref_fasta), path(ref_fai)

    output:
    tuple val(samp_name), path("*.structural_variants.vcf.gz"), val(ref_name), path("*.tbi"), emit: structural_variants

    script:
    """
    pbsv discover --hifi ${aligned_bam} ${samp_name}.svsig.gz --sample ${samp_name}
    pbsv call -j ${task.cpus} --hifi ${ref_fasta} ${samp_name}.svsig.gz ${samp_name}.${ref_name}.structural_variants.vcf
    bgzip ${samp_name}.${ref_name}.structural_variants.vcf
    tabix -p vcf ${samp_name}.${ref_name}.structural_variants.vcf.gz
    """
}

process sawfish {
    // we need to check the performance of this task for samples with <10x coverage
    label 'very_large'
    container 'quay.io/pacbio/sawfish:2.2.1_build1'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index), path(ref_fasta), path(ref_fai)

    output:
    tuple val(samp_name), path("*.structural_variants.vcf.gz"), val(ref_name), path("*.tbi"), emit: structural_variants

    script:
    """
    sawfish discover \
        --threads ${task.cpus} \
        --disable-cnv \
        --ref ${ref_fasta} \
        --bam ${aligned_bam} \
        --output-dir ${samp_name}_discover

    sawfish joint-call \
        --threads ${task.cpus} \
        --sample ${samp_name}_discover \
        --output-dir ${samp_name}_call
    
    mv -v ${samp_name}_call/genotyped.sv.vcf.gz ${samp_name}.${ref_name}.structural_variants.vcf.gz
    mv -v ${samp_name}_call/genotyped.sv.vcf.gz.tbi ${samp_name}.${ref_name}.structural_variants.vcf.gz.tbi
    """
}

process hiphase {
    publishDir "${params.outdir}/phased_output/${samp_name}", mode: 'copy'
    label 'large'
    container 'quay.io/pacbio/hiphase:1.5.0_build1'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index), path(small_variant_vcf), path(small_variant_vcf_index), path(structural_variant_vcf), path(structural_variant_vcf_index), path(ref_fasta), path(ref_fai)

    output:
    tuple val(samp_name), path("*.haplotagged.bam"), val(ref_name), path("*.haplotagged.bam.bai"), emit: hap_phased_bams
    tuple path("*.phased.vcf.gz*"), path("*.summary.tsv")

    script:
    """
    hiphase \
        --sample-name ${samp_name} \
        --threads ${task.cpus} \
        --reference ${ref_fasta} \
        --bam ${aligned_bam} \
        --output-bam ${samp_name}.${ref_name}.haplotagged.bam \
        --vcf ${small_variant_vcf} \
        --output-vcf ${samp_name}.${ref_name}.small_variants.phased.vcf.gz \
        --vcf ${structural_variant_vcf} \
        --output-vcf ${samp_name}.${ref_name}.structural_variants.phased.vcf.gz \
        --summary-file ${samp_name}.hiphase.summary.tsv
    """
}

process call_msps {
    publishDir "${params.outdir}/fire_bams/${samp_name}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index)

    output:
    tuple val(samp_name), path("*.nucs.bam"), val(ref_name), path("*.nucs.bam.bai"), emit: msp_bams

    script:
    """
    conda run -n fiberseq-qc ft add-nucleosomes \
        --threads ${task.cpus} --ml ${params.confidence_ml_val} \
        -v \
        ${aligned_bam} ${samp_name}.${ref_name}.6ma.nucs.bam;
    samtools index -@ ${task.cpus - 1} ${samp_name}.${ref_name}.6ma.nucs.bam;
    """
}

process fiberseq_qc {
    publishDir "${params.outdir}/fiberseq-qc/${samp_name}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index)

    output:
    tuple val(samp_name), path("*"), val(ref_name)

    script:
    """
    conda run -n fiberseq-qc runall-qc.V2.tcsh ./ ${samp_name} ${aligned_bam};
    awk -F '\t' '{total_A+=\$12}; {meth_A+=\$13} END {print meth_A/total_A}' *.all.tbl > ${samp_name}.6ma.methylation_rate.txt;
    """
}

process create_pileups {
    publishDir "${params.outdir}/pileups/${samp_name}", mode: 'copy'
    label 'medium'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index)

    output:
    tuple val(samp_name), path("*.tsv"), val(ref_name), emit: pileups

    script:
    """
    ft pileup \
        --m6a --ml ${params.confidence_ml_val} \
        --cpg \
        -t ${task.cpus} \
        --ftx "len(msp)>${params.minimum_msp_dist}" \
        ${aligned_bam} \
    | awk -v OFS="\t" -v FS="\t" '{print \$1,\$2,\$3,\$9/(\$4+0.1),\$10/(\$4+0.1),\$7/(\$4+0.1)}' \
    > ${samp_name}.pileup_all.tsv
    """
}

process pileupbedgraphtobigwig {
    publishDir "${params.outdir}/pileups/${samp_name}", mode: 'copy'
    label 'large'
    container 'quay.io/pacbio/bigtools:3844b58_build1'

    input:
    tuple val(samp_name), path(bedgraph), val(ref_name), path(ref_fai), val(feature), val(col_num)

    output:
    tuple val(samp_name), path("*.bw")

    script:
    """
    cut -f 1,2 ${ref_fai} > chromsizes
    cut -f 1,2,3,${col_num} ${bedgraph} | grep -v '^#' | sort -k1,1 -k2,2n > temp.bedgraph
    
    bedgraphtobigwig \
        --nthreads ${task.cpus} \
        temp.bedgraph \
        chromsizes \
        ${samp_name}.${feature}.bw
    """
}

/////////////////////////////////////////////////////////////

workflow {

    references_ch = channel.of(
            [
                "${params.ref_path}/chm13/chm13.fasta",
                "${params.ref_path}/chm13/chm13.fasta.fai",
                "chm13",
            ]
        )
        .concat(
            channel.of(
                [
                    "${params.ref_path}/hg38/hg38.fasta",
                    "${params.ref_path}/hg38/hg38.fasta.fai",
                    "hg38",
                ]
            )
        )
        .concat(
            channel.of(
                [
                    "${params.ref_path}/mm39/mm39.fasta",
                    "${params.ref_path}/mm39/mm39.fasta.fai",
                    "mm39",
                ]
            )
        )
        .concat(
            channel.of(
                [
                    "${params.ref_path}/mm10/mm10.fasta",
                    "${params.ref_path}/mm10/mm10.fasta.fai",
                    "mm10",
                ]
            )
        )
    // -> ref_fasta, ref_fai, ref_name
    if (params.debug) {
        references_ch.view { v -> "Available reference genome: ${v[2]}" }
    }

    // read in sample sheet and group by sample name
    // sample sheet columns: samp_name, bam_path, ref_name
    sample_sheet_ch = channel.fromPath("${params.sample_sheet}")
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> tuple(row[0], file(row[1]), row[2]) }
    if (params.debug) {
        sample_sheet_ch.view { v -> "Sample sheet entry: sample ${v[0]}, bam ${v[1]}, reference genome ${v[2]}" }
    }

    aligned_bams_input_ch = sample_sheet_ch
        .combine(references_ch, by: 2)
        .map { row -> tuple(row[1], row[2], row[0], row[3]) }
    // -> samp_name, bam_paths, ref_name, ref_fasta
    if (params.debug) {
        aligned_bams_input_ch.view { v -> "To be aligned: sample ${v[0]}, bam ${v[1]}, reference genome ${v[2]}" }
    }

    // align the bams to the reference genome
    align_bams(aligned_bams_input_ch)
    // -> samp_name, aligned_bam, ref_name, aligned_bam_index

    align_bams.out.aligned_bam
        .groupTuple(by: 0)
        .map { row -> tuple(row[0], row[1], row[2][0], row[3]) }
        .set { merge_bams_input_ch }
    // -> samp_name, [aligned_bam], ref_name, aligned_bam_index

    // merge the aligned bams based on sample name
    merge_bams(merge_bams_input_ch)
    // -> samp_name, merged_bam, ref_name, merged_bam_index

    if (params.pb_qc) {
        // If --pb_qc is set in command line, generate pacbio qc reports
        pacbio_qc(merge_bams.out.merged_bams)
    }
    // -> json reports, png plots

    if (params.phase_reads) {
        // If --phase_reads is set in command line, run hiphase to generate haplotype phased bams
        merge_bams.out.merged_bams
            .combine(references_ch, by: 2)
            .map { row -> tuple(row[1], row[2], row[0], row[3], row[4], row[5]) }
            .set { variantcalling_bams_input_ch }
        // -> samp_name, aligned_bam, ref_name, bam_index, ref_fasta, ref_fai

        // call snps and indels using deepvariant
        deepvariant(variantcalling_bams_input_ch)
        // -> samp_name, small_variant_vcf, ref_name, small_variant_vcf_index

        // call structural variants using sawfish
        sawfish(variantcalling_bams_input_ch)
        // -> samp_name, structural_variant_vcf, ref_name, structural_variant_vcf_index

        merge_bams.out.merged_bams
            .combine(deepvariant.out.vcfs, by: 0)
            .combine(sawfish.out.structural_variants, by: 0)
            .combine(references_ch, by: 2)
            .map { row -> tuple(row[1], row[2], row[0], row[3], row[4], row[6], row[7], row[9], row[10], row[11]) }
            .set { hiphase_input }
        // -> samp_name, aligned_bam, ref_name, bam_index, small_variant_vcf, small_variant_vcf_index, structural_variant_vcf, structural_variant_vcf_index, ref_fasta, ref_fai

        // run hiphase to generate haplotype phased bams
        hiphase(hiphase_input)
        // -> samp_name, haplotagged_bam, ref_name, haplotagged_bam_index
        call_msps_input_ch = hiphase.out.hap_phased_bams
    }
    else {
        call_msps_input_ch = merge_bams.out.merged_bams
    }

    // add nucleosomes and MSPs to the bams
    call_msps(call_msps_input_ch)
    // -> samp_name, msp_bam, ref_name, msp_bam_index

    // run stergachis fiberseq qc on bams that have been run through add-nucleosomes
    fiberseq_qc_workflow(call_msps.out.msp_bams)
    // -> samp_name, qc_files, ref_name

    if (params.create_bigwigs) {
        // If --create_bigwigs is set in command line, create pileups and bigwigs
        create_pileups(call_msps.out.msp_bams)
        // -> samp_name, pileups, ref_name
        // create a bedgraph of 6ma calling and calculate percent 6ma coverage for each base
        create_pileups.out.pileups
            .combine(references_ch, by: 2)
            .map { row -> tuple(row[1], row[2], row[0], row[4]) }
            .set { pileup_withref_ch }
        // -> samp_name, pileups, ref_name, ref_fai
        // convert the pileup bedgraph into a bigwig for downstream purposes
        pileup_withref_ch
            .map { row -> tuple(row[0], row[1], row[2], row[3], "perc6ma", 4) }
            .concat(pileup_withref_ch.map { row -> tuple(row[0], row[1], row[2], row[3], "perccpg", 5) })
            .concat(pileup_withref_ch.map { row -> tuple(row[0], row[1], row[2], row[3], "percnuc", 6) })
            .set { pileup_bedgraph_ch }
        pileupbedgraphtobigwig(pileup_bedgraph_ch)
    }
}
