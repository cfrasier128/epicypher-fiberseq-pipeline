#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_bam_path = ''
params.input_string_filter = ''
params.outdir = params.input_bam_path + '/fiberseq_output'
params.reference_genome = 'T2T'
params.confidence_ml_val = '250'
params.minimum_msp_dist = '10'
params.pb_qc = false
params.phase_reads = false
params.create_bigwigs = false
params.sample_sheet = ''

process merge_bams{
    publishDir "$params.outdir/0_Unaligned-bam/1_Merged-bams", mode: 'copy'
    cpus 8
    memory '16 GB'
    container 'cfrasier/epi-pacbio:latest'

    input:
    tuple val(input_sample_names), val(well_names), val(barcodes), val(output_sample_name), val(reference_genome), path(input_bam_files)
    output:
    tuple val(output_sample_name), path("${output_sample_name}.${barcodes.join('_')}.merged.bam"), val(reference_genome), emit: merged_bams
    script:
    """
    pbmerge ${input_bam_files.join(' ')} -o ${output_sample_name}.${barcodes.join('_')}.merged.bam
    """
}

process align_bams{
    publishDir "$params.outdir/2_Aligned-bam/1_Initial-align-temporary"
    cpus 16
    memory '32 GB'
    container 'cfrasier/epi-pacbio:latest' 

    input:
    tuple val(sampname), path(input_bam), val(reference_genome), path(ref_mmi)
    output:
    tuple val(sampname), path("*aligned.sorted.bam"), val(reference_genome), emit: aligned_bam
    script:
    """
    pbmm2 align $ref_mmi $input_bam ${sampname}.aligned.sorted.bam -j 14 --preset HIFI --sort -J 2 --log-level INFO
    """
}

process pacbio_qc{
    cpus 1
    memory '4 GB' 
    publishDir "$params.outdir/6_Seq_Stats/$sampname", mode: 'copy'
    container 'cfrasier/epi-pacbio:latest' 

    input:
    tuple val(sampname), path(raw_bam), val(reference_genome), path(bam_index)
    output:
    path("*")
    script:
    """
    pbindex $raw_bam -j $task.cpus;
    dataset create --type ConsensusReadSet --name \$(echo $sampname | cut -f 2 -d ".") ${sampname}.hifi_reads.ccsreadset.xml $raw_bam;
    runqc-reports -b --pdf-report ${sampname}.pbqc.report.pdf ${sampname}.hifi_reads.ccsreadset.xml
    """
}

process index_bams{
    publishDir "$params.outdir/2_Aligned-bam/1_Initial-align"
    cpus 1
    memory '4 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(sorted_bam), val(reference_genome)
    output:
    tuple val(sampname), path(sorted_bam), val(reference_genome), path("*.bai"), emit: bam_windex
    script:
    """
    samtools index -@ $task.cpus $sorted_bam
    """
}

process call_variants{
    publishDir "$params.outdir/2_Aligned-bam/2_Variant-calling/", mode: 'copy'
    cpus 8
    memory '16 GB'

    input:
    tuple val(sampname), path(aligned_bam), val(reference_genome), path(bam_index), path(ref_fasta), path(ref_fai)
    output:
    tuple val(sampname), path("*vcf.gz"), val(reference_genome), path("*vcf.gz.tbi"), emit:vcfs
    script:
    """
    docker run -v \$PWD:\$PWD \
    -v /media/:/media/ \
    -w \$PWD \
    google/deepvariant:1.8.0  \
    /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO \
    --ref=$ref_fasta \
    --reads $aligned_bam \
    --output_vcf=${sampname}.vcf.gz \
    --sample_name $sampname \
    --num_shards=$task.cpus
    """
}

process call_structural_variants{
    publishDir "$params.outdir/2_Aligned-bam/2_Variant-calling/", mode: 'copy'
    cpus 8
    memory '16 GB'
    container 'quay.io/pacbio/pbsv:2.11.0_build1'

    input:
    tuple val(sampname), path(aligned_bam), val(reference_genome), path(bam_index), path(ref_fasta), path(ref_fai)
    output:
    tuple val(sampname), path("*.sv.vcf.gz"), val(reference_genome), path("*.tbi"), emit: structural_variants
    script:
    """
    pbsv discover --hifi $aligned_bam ${sampname}.svsig.gz --sample $sampname ;
    pbsv call -j $task.cpus --hifi $ref_fasta ${sampname}.svsig.gz ${sampname}.sv.vcf;
    bgzip ${sampname}.sv.vcf;
    tabix -p vcf ${sampname}.sv.vcf.gz;
    """
}

process hiphase{
    publishDir "$params.outdir/2_Aligned-bam/3_Haplotype-phased/", mode: 'copy'
    cpus 8
    memory '16 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(aligned_bam), val(reference_genome), path(bam_index), path(vcf), path(vcf_index), path(struct_vars), path(struct_tbi), path(ref_fasta), path(ref_fai)
    output:
    tuple val(sampname), path("*.haplotagged.bam"), val(reference_genome), path("*.haplotagged.bam.bai"), emit: hap_phased_bams
    path("*")
    script:
    """
    conda run -n fiberseq-qc hiphase --reference $ref_fasta \
    --threads $task.cpus --disable-global-realignment \
    --bam $aligned_bam \
    --output-bam ${sampname}.haplotagged.bam \
    --vcf $vcf \
    --output-vcf ${sampname}.haplotagged.vcf \
    --vcf $struct_vars \
    --output-vcf ${sampname}.haplotagged.structural.vcf \
    --ignore-read-groups \
    --summary-file ${sampname}.summary.tsv
    """
}

process call_msps {
    publishDir "$params.outdir/3_Fibertools/1_FIRE-bams/", mode: 'copy'
    cpus 8
    memory '16 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(aligned_bam), val(reference_genome), path(bam_index)
    output:
    tuple val(sampname), path("*.sorted.bam"), val(reference_genome), emit: msp_bams
    path("*.bai"), emit: msp_bams_index
    script:
    """
    conda run -n fiberseq-qc ft add-nucleosomes \
    --threads $task.cpus --ml $params.confidence_ml_val \
    -v \
    $aligned_bam ${sampname}.6ma.nucs.sorted.bam;
    samtools index -@ $task.cpus ${sampname}.6ma.nucs.sorted.bam;
    """
}

process fiberseq_qc {
    cpus 8
    memory '16 GB' 
    publishDir "$params.outdir/3_Fibertools/2_Fiberseq-qc/", mode: 'copy'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(aligned_bam), val(reference_genome)
    output:
    tuple val(sampname), path("*"), val(reference_genome)
    script:
    """
    conda run -n fiberseq-qc runall-qc.V2.tcsh ./ $sampname $aligned_bam;
    awk -F '\t' '{total_A+=\$12}; {meth_A+=\$13} END {print meth_A/total_A}' *.all.tbl > ${sampname}.6ma.methylation_rate.txt;
    """
}

process create_pileups {
    publishDir "$params.outdir/4_Pileups_Bigwigs/1_Pileups/"
    cpus 4
    memory '8 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(aligned_bam), val(reference_genome)
    path(bam_index)
    output:
    tuple val(sampname), path("*.tsv"), val(reference_genome), emit: pileups
    script:
    """
    ft pileup \
    --m6a --ml $params.confidence_ml_val \
    --cpg \
    -t $task.cpus \
    --ftx "len(msp)>$params.minimum_msp_dist" \
    $aligned_bam | awk -v OFS="\t" -v FS="\t" '{print \$1,\$2,\$3,\$9/(\$4+0.1),\$10/(\$4+0.1),\$7/(\$4+0.1)}' > ${sampname}.pileup_all.tsv
    """
}

process pileupbedgraphtobigwig_6ma{
    publishDir "$params.outdir/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(pileup), val(reference_genome), path(genome_chromsizes)
    output:
    tuple val(sampname), path("*.bw")
    script:
    """
    cut -f 1,2,3,4 $pileup > temp.bedgraph
    bedGraphToBigWig temp.bedgraph $genome_chromsizes ${sampname}.perc6ma.bw
    """
}

process pileupbedgraphtobigwig_5mC{
    publishDir "$params.outdir/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(pileup), val(reference_genome), path(genome_chromsizes)
    output:
    tuple val(sampname), path("*.bw")
    script:
    """
    cut -f 1-3,5 $pileup > temp.bedgraph
    bedGraphToBigWig temp.bedgraph $genome_chromsizes ${sampname}.perc5mc.bw
    """
}

process pileupbedgraphtobigwig_nuc{
    publishDir "$params.outdir/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4 GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(pileup), val(reference_genome), path(genome_chromsizes)
    output:
    tuple val(sampname), path("*.bw")
    script:
    """
    cut -f 1-3,6 $pileup > temp.bedgraph
    bedGraphToBigWig temp.bedgraph $genome_chromsizes ${sampname}.percnuc.bw
    """
}

/////////////////////////////////////////////////////////////

workflow  {
    // Set the index files based on the reference genome provided

    T2T_ref_ch = channel.of(["T2T",
                    "/media/genomics/18Tb_1/references/T2T/T2T.ch13v2.0.mmi",
                    "/media/genomics/18Tb_1/references/T2T/hs1.chrom.sizes",
                    "/media/genomics/18Tb_1/references/T2T/chm13v2.0.clean.fasta",
                    "/media/genomics/18Tb_1/references/T2T/chm13v2.0.clean.fasta.fai"])
    hg38_ref_ch = channel.of(["hg38",
                    "/media/genomics/18Tb_1/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.mmi",
                    "/media/genomics/18Tb_1/references/hg38/hg38.uscsnames.v40.NO_ALTS.chrom.sizes",
                    "/media/genomics/18Tb_1/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.fa",
                    "/media/genomics/18Tb_1/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.fa.fai"])
    mm39_ref_ch = channel.of(["mm39",
                    "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.mmi",
                    "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.chrom.sizes",
                    "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.fa",
                    "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.fa.fai"])
    
    T2T_ref_ch.concat(hg38_ref_ch).concat(mm39_ref_ch).map{row ->
        tuple( row[1], row[2], row[0], row[3], row[4] )
    }
        .set { references_ch }
    // references_ch.view()

    called_input_ch = channel.fromPath("${params.input_bam_path}/*${params.input_string_filter}*.bam")       // grabs all bams in the input path
    called_input_ch.map { file -> tuple(file.baseName.split('\\.')[1], file) }.set { input_bams_names_ch }   // get sample name, i.e. drop all extensions

    samplesheet_ch = channel.fromPath("${params.sample_sheet}")                                             // read in sample sheet
        .splitCsv(skip:1, sep: '\t')

    samplesheet_ch.combine(input_bams_names_ch, by: 0)                                                      // combine sample sheet with input bams based on sample name 
        .set { combined_input_samplesheet_ch }

    combined_input_samplesheet_ch.branch { row ->                                                           // branch based on whether to merge bams or not                  
        merge: row[3] != ""
        no_merge: row[3] == ""
    }.set { branched_samplesheet_ch }

    branched_samplesheet_ch.merge                                                                           // for samples to be merged, group based on sample name
        .groupTuple(by: 3)
        .set { merged_samples_ch }

    merge_bams(merged_samples_ch)                                                                           // Merge bams based on sample sheet info

    pre_col_ch = merge_bams.out.merged_bams                             

    pre_col_ch.map { row -> tuple( row[0], row[1], row[2][0] ) }                                            // flatten out reference genome after merge
        .set { final_merged_bams_ch }

    branched_samplesheet_ch.no_merge                                                                        // set unmerged channel to correct format to match merged bams channel
        .map { row -> tuple( row[0], row[5], row[4] ) }
        .set { no_merge_bams_ch }
    
    final_merged_bams_ch.concat(no_merge_bams_ch).set { all_bams_ch }                                       // combine merged and unmerged bams into single channel

    all_bams_ch.combine(references_ch, by: 2)                                                               // combine bams with reference genome info based on reference genome name
        .map { row -> tuple( row[1], row[2], row[0], row[3] ) }
        .set { aligned_bams_input_ch }

    align_bams(aligned_bams_input_ch)                                                                       // align bams
    index_bams(align_bams.out.aligned_bam)

    if (params.pb_qc){                                                                                      // If --pb_qc is set in command line, generate pacbio qc reports
        pacbio_qc(index_bams.out.bam_windex)                                                              
    }

    if (params.phase_reads){                                                                                // If --phase_reads is set in command line, run hiphase to generate haplotype phased bams
        index_bams.out.bam_windex.combine(references_ch, by: 2)                                             // combine bams with reference genome info based on reference genome name
        .map { row -> tuple( row[1], row[2], row[0], row[3], row[6], row[7] ) }
        .set { variantcalling_bams_input_ch }
        call_variants (variantcalling_bams_input_ch)                                                        // call snps and indels using deepvariant
        call_structural_variants(variantcalling_bams_input_ch)                                              // call structural variants using pbsv
        index_bams.out.bam_windex.combine(call_variants.out.vcfs, by: 0).set{bams_snps}                     // combine the bams and vcf channels
        bams_snps.combine(call_structural_variants.out.structural_variants, by: 0).set{hiphase_input}       // combine the bams and sv channels
        hiphase_input.combine(references_ch, by: 2)                                             // combine bams with reference genome info based on reference genome name
        .map { row -> tuple( row[1], row[2], row[0], row[3], row[4], row[6], row[7], row[9], row[12], row[13] ) }
        .set { hiphase_input_withref }
        hiphase(hiphase_input_withref)
        call_msps_input_ch = hiphase.out.hap_phased_bams
    }
    else {
        call_msps_input_ch = index_bams.out.bam_windex
    }

    call_msps(call_msps_input_ch)                                                                            // add nucleosomes and MSPs to the bams
    fiberseq_qc(call_msps.out.msp_bams)                                                                      // run stergachis fiberseq qc on bams that have been run through add-nucleosomes


    if (params.create_bigwigs){                                                                              // If --create_bigwigs is set in command line, create pileups and bigwigs
        create_pileups(call_msps.out.msp_bams,call_msps.out.msp_bams_index)                                  // create a bedgraph of 6ma calling and calculate percent 6ma coverage for each base
        create_pileups.out.pileups.combine(references_ch, by: 2)                                       // combine pileup with reference genome info based on reference genome name
        .map { row -> tuple( row[1], row[2], row[0], row[4] ) }
        .set { pileup_withref_ch }
        pileup_withref_ch.view()  
        pileupbedgraphtobigwig_6ma(pileup_withref_ch)                                // convert the pileup bedgraph into a bigwig for downstream purposes
        pileupbedgraphtobigwig_5mC(pileup_withref_ch)                                // convert the pileup bedgraph into a bigwig for downstream purposes
        pileupbedgraphtobigwig_nuc(pileup_withref_ch)                                // convert the pileup bedgraph into a bigwig for downstream purposes
    }
}