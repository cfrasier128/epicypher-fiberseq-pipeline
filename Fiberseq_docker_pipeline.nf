#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// inputs
params.sample_sheet = ''
params.ref_path = ''

// default parameters
params.ref_name = 'T2T'
params.confidence_ml_val = '250'
params.minimum_msp_dist = '10'

// optional steps
params.pb_qc = false
params.phase_reads = false
params.create_bigwigs = false

// output directory
params.outdir = "${workflow.launchDir}/fiberseq_output"

process merge_bams {
    publishDir "${params.outdir}/0_Unaligned-bam/1_Merged-bams", mode: 'copy'
    cpus 8
    memory '16 GB'
    container 'cfrasier/epi-pacbio:latest'

    input:
    tuple val(input_sample_names), val(well_names), val(barcodes), val(output_sample_name), val(ref_name), path(input_bam_files)

    output:
    tuple val(output_sample_name), path("${output_sample_name}.${barcodes.join('_')}.merged.bam"), val(ref_name), emit: merged_bams

    script:
    """
    pbmerge ${input_bam_files.join(' ')} -o ${output_sample_name}.${barcodes.join('_')}.merged.bam
    """
}

process align_bams {
    publishDir "${params.outdir}/2_Aligned-bam/1_Initial-align-temporary"
    cpus 16
    memory '32 GB'
    container 'quay.io/pacbio/pbmm2:1.17.0_build1'

    input:
    tuple val(samp_name), path(input_bam), val(ref_name), path(ref_fasta)

    output:
    tuple val(samp_name), path("*${ref_name}.bam"), val(ref_name), path("*${ref_name}.bam.bai"), emit: aligned_bam

    script:
    """
    pbmm2 align \
        --preset HIFI --log-level INFO \
        --num-threads ${task.cpus} \
        --sort --sort-memory 4G --bam-index BAI \
        --sample ${samp_name} \
        ${ref_fasta} \
        ${input_bam} \
        ${samp_name}.${ref_name}.bam
    """
}

process pacbio_qc {
    cpus 1
    memory '4 GB'
    publishDir "${params.outdir}/6_Seq_Stats/${samp_name}", mode: 'copy'
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
    publishDir "${params.outdir}/2_Aligned-bam/2_Variant-calling/", mode: 'copy'
    cpus 16
    memory '64 GB'
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
    publishDir "${params.outdir}/2_Aligned-bam/2_Variant-calling/", mode: 'copy'
    cpus 8
    memory '16 GB'
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
    publishDir "${params.outdir}/2_Aligned-bam/2_Variant-calling/", mode: 'copy'
    cpus 8
    memory '64 GB'
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
    publishDir "${params.outdir}/2_Aligned-bam/3_Haplotype-phased/", mode: 'copy'
    cpus 8
    memory '16 GB'
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
    publishDir "${params.outdir}/3_Fibertools/1_FIRE-bams/", mode: 'copy'
    cpus 8
    memory '16 GB'
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
    cpus 8
    memory '16 GB'
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
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
    publishDir "${params.outdir}/4_Pileups_Bigwigs/1_Pileups/"
    cpus 4
    memory '8 GB'
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
    publishDir "${params.outdir}/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(pileup), val(ref_name), path(ref_fai), val(feature), val(col_num)

    output:
    tuple val(samp_name), path("*.bw")

    script:
    """
    cut -f 1,2 ${ref_fai} > chromsizes
    cut -f 1,2,3,${col_num} ${pileup} > temp.bedgraph
    bedGraphToBigWig temp.bedgraph chromsizes ${samp_name}.${feature}.bw
    """
}

/////////////////////////////////////////////////////////////

workflow {

    references_ch = channel.of(
            [
                "${params.ref_path}/T2T/chm13v2.0.clean.fasta",
                "${params.ref_path}/T2T/chm13v2.0.clean.fasta.fai",
                "T2T",
            ]
        )
        .concat(
            channel.of(
                [
                    "${params.ref_path}/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.fa",
                    "${params.ref_path}/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.fa.fai",
                    "hg38",
                ]
            )
        )
        .concat(
            channel.of(
                [
                    "${params.ref_path}/mm10/GCF_000001635.27_GRCm39_genomic.fa",
                    "${params.ref_path}/mm10/GCF_000001635.27_GRCm39_genomic.fa.fai",
                    "mm39",
                ]
            )
        )
    // -> ref_fasta, ref_fai, ref_name
    // references_ch.view { v -> "Available reference genome: ${v[2]}" }

    // read in sample sheet and group by sample name
    // sample sheet columns: samp_name, bam_path, ref_name
    samplesheet_ch = channel.fromPath("${params.sample_sheet}")
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> tuple(row[0], file(row[1]), row[2]) }
        .groupTuple(by: 0)
        .branch { row ->
            merge: row[1].size() > 1
            no_merge: row[1].size() == 1
        }
    // -> samp_name, [bam_paths], [ref_name]
    // branch based on whether to merge or not based on if samplename is provided in sample sheet
    // samplesheet_ch.view{ v -> "Sample from sheet: ${v[0]}, bam(s): ${v[1].join(', ')}, reference genome: ${v[2]}" }
    // samplesheet_ch.merge.view { v -> "For sample ${v[0]}, merging BAMs ${v[1].join(', ')}" }
    // samplesheet_ch.no_merge.view { v -> "For sample ${v[0]}, no merging needed for BAM ${v[1][0]}" }

    // merge bams based on sample sheet info
    merge_bams(samplesheet_ch.merge).merged_bams.map { row -> tuple(row[0], row[1], row[2][0]) }.set { final_merged_bams_ch }
    // -> samp_name, merged_bam, ref_name

    // for samples that do not need merging, just pass through the bam paths
    no_merge_bams_ch = samplesheet_ch.no_merge.map { row -> tuple(row[0], row[1][0], row[2][0]) }
    // -> samp_name, bam_path, ref_name

    // combine merged and unmerged bams into single channel
    all_bams_ch = final_merged_bams_ch.concat(no_merge_bams_ch)
    // -> samp_name, bam_path, ref_name
    // all_bams_ch.view { v -> "To be aligned: sample ${v[0]}, bam ${v[1]}, reference genome ${v[2]}" }

    // align the bams to the reference genome
    all_bams_ch
        .combine(references_ch, by: 2)
        .map { row -> tuple(row[1], row[2], row[0], row[3]) }
        .set { aligned_bams_input_ch }
    // -> samp_name, bam_path, ref_name, ref_fasta
    align_bams(aligned_bams_input_ch)
    // -> samp_name, aligned_bam, ref_name, aligned_bam_index

    if (params.pb_qc) {
        // If --pb_qc is set in command line, generate pacbio qc reports
        pacbio_qc(align_bams.out.aligned_bam)
    }
    // -> json reports, png plots

    if (params.phase_reads) {
        // If --phase_reads is set in command line, run hiphase to generate haplotype phased bams
        align_bams.out.aligned_bam
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

        align_bams.out.aligned_bam
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
        call_msps_input_ch = align_bams.out.aligned_bam
    }

    // add nucleosomes and MSPs to the bams
    call_msps(call_msps_input_ch)
    // -> samp_name, msp_bam, ref_name, msp_bam_index

    // run stergachis fiberseq qc on bams that have been run through add-nucleosomes
    fiberseq_qc(call_msps.out.msp_bams)
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
        pileup_withref_ch.map { row -> tuple(row[0], row[1], row[2], row[3], "perc6ma", 4) }
            .concat(pileup_withref_ch.map { row -> tuple(row[0], row[1], row[2], row[3], "perccpg", 5) })
            .concat(pileup_withref_ch.map { row -> tuple(row[0], row[1], row[2], row[3], "percnuc", 6) })
            .set { pileup_bedgraph_ch }
        pileupbedgraphtobigwig(pileup_bedgraph_ch)
    }
}
