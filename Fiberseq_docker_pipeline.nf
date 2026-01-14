#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// inputs
params.input_bam_path = ''
params.input_string_filter = ''
params.sample_sheet = ''
params.ref_path = '/media/genomics/18Tb_1/references'

// default parameters
params.ref_name = 'T2T'
params.confidence_ml_val = '250'
params.minimum_msp_dist = '10'

// optional steps
params.pb_qc = false
params.phase_reads = false
params.create_bigwigs = false

// output directory
params.outdir = params.input_bam_path + '/fiberseq_output'

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
    container 'cfrasier/epi-pacbio:latest'

    input:
    tuple val(samp_name), path(input_bam), val(ref_name), path(ref_fasta)

    output:
    tuple val(samp_name), path("*aligned.sorted.bam"), val(ref_name), emit: aligned_bam

    script:
    """
    pbmm2 align ${ref_fasta} ${input_bam} ${samp_name}.aligned.sorted.bam -j 14 --preset HIFI --sort -J 2 --log-level INFO
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
    path("ccs.report.json"), path("*plot.png"), emit: pacbio_qc_reports

    script:
    """
    pbindex ${raw_bam} -j ${task.cpus};
    dataset create --type ConsensusReadSet --name \$(echo ${samp_name} | cut -f 2 -d ".") ${samp_name}.hifi_reads.ccsreadset.xml ${raw_bam};
    runqc-reports -b --pdf-report ${samp_name}.pbqc.report.pdf ${samp_name}.hifi_reads.ccsreadset.xml
    """
}

process index_bams {
    publishDir "${params.outdir}/2_Aligned-bam/1_Initial-align"
    cpus 1
    memory '4 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(sorted_bam), val(ref_name)

    output:
    tuple val(samp_name), path(sorted_bam), val(ref_name), path("*.bai"), emit: bam_windex

    script:
    """
    samtools index -@ ${task.cpus} ${sorted_bam}
    """
}

process call_variants {
    publishDir "${params.outdir}/2_Aligned-bam/2_Variant-calling/", mode: 'copy'
    cpus 8
    memory '16 GB'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index), path(ref_fasta), path(ref_fai)

    output:
    tuple val(samp_name), path("*vcf.gz"), val(ref_name), path("*vcf.gz.tbi"), emit: vcfs

    script:
    """
    docker run -v \$PWD:\$PWD \
    -v /media/:/media/ \
    -w \$PWD \
    google/deepvariant:1.8.0  \
    /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO \
        --ref=${ref_fasta} \
        --reads ${aligned_bam} \
        --output_vcf=${samp_name}.vcf.gz \
        --sample_name ${samp_name} \
        --num_shards=${task.cpus}
    """
}

process call_structural_variants {
    publishDir "${params.outdir}/2_Aligned-bam/2_Variant-calling/", mode: 'copy'
    cpus 8
    memory '16 GB'
    container 'quay.io/pacbio/pbsv:2.11.0_build1'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index), path(ref_fasta), path(ref_fai)

    output:
    tuple val(samp_name), path("*.sv.vcf.gz"), val(ref_name), path("*.tbi"), emit: structural_variants

    script:
    """
    pbsv discover --hifi ${aligned_bam} ${samp_name}.svsig.gz --sample ${samp_name} ;
    pbsv call -j ${task.cpus} --hifi ${ref_fasta} ${samp_name}.svsig.gz ${samp_name}.sv.vcf;
    bgzip ${samp_name}.sv.vcf;
    tabix -p vcf ${samp_name}.sv.vcf.gz;
    """
}

process hiphase {
    publishDir "${params.outdir}/2_Aligned-bam/3_Haplotype-phased/", mode: 'copy'
    cpus 8
    memory '16 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name), path(bam_index), path(vcf), path(vcf_index), path(struct_vars), path(struct_tbi), path(ref_fasta), path(ref_fai)

    output:
    tuple val(samp_name), path("*.haplotagged.bam"), val(ref_name), path("*.haplotagged.bam.bai"), emit: hap_phased_bams
    path "*"

    script:
    """
    conda run -n fiberseq-qc hiphase --reference ${ref_fasta} \
        --threads ${task.cpus} --disable-global-realignment \
        --bam ${aligned_bam} \
        --output-bam ${samp_name}.haplotagged.bam \
        --vcf ${vcf} \
        --output-vcf ${samp_name}.haplotagged.vcf \
        --vcf ${struct_vars} \
        --output-vcf ${samp_name}.haplotagged.structural.vcf \
        --ignore-read-groups \
        --summary-file ${samp_name}.summary.tsv
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
    tuple val(samp_name), path("*.sorted.bam"), val(ref_name), emit: msp_bams
    path ("*.bai"), emit: msp_bams_index

    script:
    """
    conda run -n fiberseq-qc ft add-nucleosomes \
        --threads ${task.cpus} --ml ${params.confidence_ml_val} \
        -v \
        ${aligned_bam} ${samp_name}.6ma.nucs.sorted.bam;
    samtools index -@ ${task.cpus} ${samp_name}.6ma.nucs.sorted.bam;
    """
}

process fiberseq_qc {
    cpus 8
    memory '16 GB'
    publishDir "${params.outdir}/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(aligned_bam), val(ref_name)

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
    tuple val(samp_name), path(aligned_bam), val(ref_name)
    path bam_index

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

process pileupbedgraphtobigwig_6ma {
    publishDir "${params.outdir}/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(pileup), val(ref_name), path(ref_fai)

    output:
    tuple val(samp_name), path("*.bw")

    script:
    """
    cut -f 1,2 ${ref_fai} > chromsizes
    cut -f 1,2,3,4 ${pileup} > temp.bedgraph
    bedGraphToBigWig temp.bedgraph chromsizes ${samp_name}.perc6ma.bw
    """
}

process pileupbedgraphtobigwig_5mC {
    publishDir "${params.outdir}/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(pileup), val(ref_name), path(ref_fai)

    output:
    tuple val(samp_name), path("*.bw")

    script:
    """
    cut -f 1,2 ${ref_fai} > chromsizes
    cut -f 1-3,5 ${pileup} > temp.bedgraph
    bedGraphToBigWig temp.bedgraph chromsizes ${samp_name}.perc5mc.bw
    """
}

process pileupbedgraphtobigwig_nuc {
    publishDir "${params.outdir}/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(samp_name), path(pileup), val(ref_name), path(ref_fai)

    output:
    tuple val(samp_name), path("*.bw")

    script:
    """
    cut -f 1,2 ${ref_fai} > chromsizes
    cut -f 1-3,6 ${pileup} > temp.bedgraph
    bedGraphToBigWig temp.bedgraph chromsizes ${samp_name}.percnuc.bw
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
    .concat(channel.of(
        [
            "${params.ref_path}/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.fa",
            "${params.ref_path}/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.fa.fai",
            "hg38",
        ]
    ))
    .concat(channel.of(
        [
            "${params.ref_path}/mm10/GCF_000001635.27_GRCm39_genomic.fa",
            "${params.ref_path}/mm10/GCF_000001635.27_GRCm39_genomic.fa.fai",
            "mm39",
        ]
    ))
    // -> ref_fasta, ref_fai, ref_name

    called_input_ch = channel.fromPath("${params.input_bam_path}/*${params.input_string_filter}*.bam")
    // grabs all bams in the input path
    called_input_ch.map { file -> tuple(file.baseName.split('\\.')[1], file) }.set { input_bams_names_ch }
    // get sample name, i.e. drop all extensions; if file is {movie}.{sample_name}.bam, get sample_name, assumes no other dots in sample name
    // -> sample_name, path_to_bam

    samplesheet_ch = channel.fromPath("${params.sample_sheet}")
        .splitCsv(skip: 1, sep: '\t')
    // -> biosamplename, wellname, barcode, samplename, ref_name, path

    samplesheet_ch
        .combine(input_bams_names_ch, by: 0)
        .set { combined_input_samplesheet_ch }
    // -> biosamplename, wellname, barcode, samplename, ref_name, path, sample_name, path_to_bam

    combined_input_samplesheet_ch
        .branch { row ->
            merge: row[3] != ""
            no_merge: row[3] == ""
        }
        .set { branched_samplesheet_ch }
    // branch based on whether to merge or not based on if samplename is provided in sample sheet

    branched_samplesheet_ch.merge
        .groupTuple(by: 3)
        .set { merged_samples_ch }
    // -> biosamplename, wellname, barcode, samplename, ref_name, path
    // group by samplename to prepare for merging bams

    merge_bams(merged_samples_ch)
    // Merge bams based on sample sheet info
    // -> samp_name, merged_bam, ref_name

    merge_bams.out.merged_bams
        .map { row -> tuple(row[0], row[1], row[2][0]) }
        .set { final_merged_bams_ch }
    // -> samp_name, merged_bam, ref_name
    // take first reference genome from list

    branched_samplesheet_ch.no_merge
        .map { row -> tuple(row[0], row[5], row[4]) }
        .set { no_merge_bams_ch }
    // -> samp_name, bam_path, ref_name

    final_merged_bams_ch.concat(no_merge_bams_ch).set { all_bams_ch }
    // combine merged and unmerged bams into single channel
    // -> samp_name, bam_path, ref_name

    all_bams_ch
        .combine(references_ch, by: 2)
        .map { row -> tuple(row[1], row[2], row[0], row[3]) }
        .set { aligned_bams_input_ch }
    // -> samp_name, bam_path, ref_name, ref_fasta

    align_bams(aligned_bams_input_ch)
    // -> samp_name, aligned_bam, ref_name

    index_bams(align_bams.out.aligned_bam)
    // -> samp_name, aligned_bam, ref_name, bam_index

    if (params.pb_qc) {
        // If --pb_qc is set in command line, generate pacbio qc reports
        pacbio_qc(index_bams.out.bam_windex)
    }
    // -> json reports, png plots

    if (params.phase_reads) {
        // If --phase_reads is set in command line, run hiphase to generate haplotype phased bams
        index_bams.out.bam_windex
            .combine(references_ch, by: 2)
            .map { row -> tuple(row[1], row[2], row[0], row[3], row[4], row[5]) }
            .set { variantcalling_bams_input_ch }
        // -> samp_name, aligned_bam, ref_name, bam_index, ref_fasta, ref_fai

        // call snps and indels using deepvariant
        call_variants(variantcalling_bams_input_ch)
        // -> samp_name, small_variant_vcf, ref_name, small_variant_vcf_index

        // call structural variants using pbsv
        call_structural_variants(variantcalling_bams_input_ch)
        // -> samp_name, structural_variant_vcf, ref_name, structural_variant_vcf_index

        index_bams.out.bam_windex.combine(call_variants.out.vcfs, by: 0).set { bams_snps }
        // combine the bams and vcf channels
        // -> samp_name, aligned_bam, ref_name, bam_index, small_variant_vcf, small_variant_vcf_index

        bams_snps.combine(call_structural_variants.out.structural_variants, by: 0).set { hiphase_input }
        // combine the bams and sv channels
        // -> samp_name, aligned_bam, ref_name, bam_index, small_variant_vcf, small_variant_vcf_index, structural_variant_vcf, structural_variant_vcf_index

        hiphase_input
            .combine(references_ch, by: 2)
            .map { row -> tuple(row[1], row[2], row[0], row[3], row[4], row[6], row[7], row[9], row[10], row[11]) }
            .set { hiphase_input_withref }
        // -> samp_name, aligned_bam, ref_name, bam_index, small_variant_vcf, small_variant_vcf_index, structural_variant_vcf, structural_variant_vcf_index, ref_fasta, ref_fai
        
        hiphase(hiphase_input_withref)
        // -> samp_name, haplotagged_bam, ref_name, haplotagged_bam_index
        call_msps_input_ch = hiphase.out.hap_phased_bams
    }
    else {
        call_msps_input_ch = index_bams.out.bam_windex
    }

    call_msps(call_msps_input_ch)
    // add nucleosomes and MSPs to the bams
    // -> samp_name, msp_bam, ref_name
    // -> msp_bam_index

    fiberseq_qc(call_msps.out.msp_bams)
    // run stergachis fiberseq qc on bams that have been run through add-nucleosomes
    // -> samp_name, qc_files, ref_name


    if (params.create_bigwigs) {
        // If --create_bigwigs is set in command line, create pileups and bigwigs
        create_pileups(call_msps.out.msp_bams, call_msps.out.msp_bams_index)
        // -> samp_name, pileups, ref_name
        // create a bedgraph of 6ma calling and calculate percent 6ma coverage for each base
        create_pileups.out.pileups
            .combine(references_ch, by: 2)
            .map { row -> tuple(row[1], row[2], row[0], row[4]) }
            .set { pileup_withref_ch }
        // -> samp_name, pileups, ref_name, ref_fai
        //pileup_withref_ch.view()
        // convert the pileup bedgraph into a bigwig for downstream purposes
        pileupbedgraphtobigwig_6ma(pileup_withref_ch)
        pileupbedgraphtobigwig_5mC(pileup_withref_ch)
        pileupbedgraphtobigwig_nuc(pileup_withref_ch)
        // -> samp_name, bigwig
    }
}
