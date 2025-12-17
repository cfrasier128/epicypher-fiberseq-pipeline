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

process align_bams{
    publishDir "$params.outdir/2_Aligned-bam/1_Initial-align-temporary"
    cpus 16
    memory '32 GB'
    container 'cfrasier/epi-pacbio:latest' 

    input:
    tuple val(sampname), path(input_bam)
    path(ref_mmi)
    output:
    tuple val(sampname), path("*aligned.sorted.bam"), emit: aligned_bam
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
    tuple val(sampname), path(raw_bam), path(bam_index)
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
    tuple val(sampname), path(sorted_bam)
    output:
    tuple val(sampname), path(sorted_bam), path("*.bai"), emit: bam_windex
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
    tuple val(sampname), path(aligned_bam), path(bam_index)
    path(ref_fasta)
    path(ref_fai)
    output:
    tuple val(sampname), path("*vcf.gz"), path("*vcf.gz.tbi"), emit:vcfs
    script:
    """
    docker run -v \$PWD:\$PWD \
    -v /media/:/media/ \
    -w \$PWD \
    google/deepvariant:1.8.0  \
    /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO \
    --regions "NC_000001.11 NC_000002.12 NC_000003.12 NC_000004.12 NC_000005.10 NC_000006.12 NC_000007.14 NC_000008.11 NC_000009.12 NC_000010.11 NC_000011.10 NC_000012.12 NC_000013.11 NC_000014.9 NC_000015.10 NC_000016.10 NC_000017.11 NC_000018.10 NC_000019.10 NC_000020.11 NC_000021.9 NC_000022.11 NC_000023.11 NC_000024.10" \
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
    tuple val(sampname), path(aligned_bam), path(bam_index)
    path(ref_fasta)
    path(ref_fai)
    output:
    tuple val(sampname), path("*.sv.vcf.gz"), path("*.tbi"), emit: structural_variants
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
    tuple val(sampname), path(aligned_bam), path(bam_index), path(vcf), path(vcf_index), path(struct_vars), path(tbi)
    path(ref_fasta)
    path(ref_fai)
    output:
    tuple val(sampname), path("*.haplotagged.bam"), path("*.haplotagged.bam.bai"), emit: hap_phased_bams
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
    tuple val(sampname), path(aligned_bam), path(bam_index)
    output:
    tuple val(sampname), path("*.sorted.bam"), emit: msp_bams
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
    container '`'

    input:
    tuple val(sampname), path(aligned_bam)
    output:
    tuple val(sampname), path("*")
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
    tuple val(sampname), path(aligned_bam)
    path(bam_index)
    output:
    tuple val(sampname), path("*.tsv"), emit: pileups
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
    tuple val(sampname), path(pileup)
    path(genome_chromsizes)
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
    tuple val(sampname), path(pileup)
    path(genome_chromsizes)
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
    tuple val(sampname), path(pileup)
    path(genome_chromsizes)
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
    if (params.reference_genome == 'T2T') {
        ref_mmi = "/media/genomics/18Tb_1/references/T2T/T2T.ch13v2.0.mmi"
        ref_chromsizes = "/media/genomics/18Tb_1/references/T2T/hs1.chrom.sizes"
        ref_fasta = "/media/genomics/18Tb_1/references/T2T/chm13v2.0.clean.fasta"
        ref_fai = "/media/genomics/18Tb_1/references/T2T/chm13v2.0.clean.fasta.fai"
    }
    else if (params.reference_genome == 'hg38') {
        ref_mmi = "/media/genomics/18Tb_1/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.mmi"
        ref_chromsizes = "/media/genomics/18Tb_1/references/hg38/hg38.uscsnames.v40.NO_ALTS.chrom.sizes"
        ref_fasta = "/media/genomics/18Tb_1/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.fa"
        ref_fai = "/media/genomics/18Tb_1/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.NO_ALTS.fa.fai"
    }
    else if (params.reference_genome == 'mm39') {
        ref_mmi = "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.mmi"
        ref_chromsizes = "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.chrom.sizes"
        ref_fasta = "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.fa"
        ref_fai = "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.fa.fai"
    }
    else {
        exit 1
    }

    called_input_ch = channel.fromPath("${params.input_bam_path}/*${params.input_string_filter}*.bam")       // grabs all bams in the input path
    called_input_ch.map { file -> tuple(file.baseName.split('\\.')[0], file) }.set { input_bams_names_ch }   // get sample name, i.e. drop all extensions
    input_bams_names_ch.view()
                                                                                                    
    align_bams(input_bams_names_ch,ref_mmi)                                                                  // align bams
    index_bams(align_bams.out.aligned_bam)

    if (params.pb_qc){                                                                                       // If --pb_qc is set in command line, generate pacbio qc reports
        pacbio_qc(index_bams.out.bam_windex)                                                              
    }

    if (params.phase_reads){                                                                                 // If --phase_reads is set in command line, run hiphase to generate haplotype phased bams
        call_variants (index_bams.out.bam_windex,ref_fasta,ref_fai)                                          // call snps and indels using deepvariant
        call_structural_variants(index_bams.out.bam_windex,ref_fasta,ref_fai)                                // call structural variants using pbsv
        index_bams.out.bam_windex.combine(call_variants.out.vcfs, by: 0).set{bams_snps}                      // combine the bams and vcf channels
        bams_snps.combine(call_structural_variants.out.structural_variants, by: 0).set{hiphase_input}        // combine the bams and sv channels
        hiphase(hiphase_input,ref_fasta,ref_fai)
        call_msps_input_ch = hiphase.out.hap_phased_bams
    }
    else {
        call_msps_input_ch = index_bams.out.bam_windex
    }

    call_msps(call_msps_input_ch)                                                                            // add nucleosomes and MSPs to the bams
    fiberseq_qc(call_msps.out.msp_bams)                                                                      // run stergachis fiberseq qc on bams that have been run through add-nucleosomes


    if (params.create_bigwigs){                                                                              // If --create_bigwigs is set in command line, create pileups and bigwigs
        create_pileups(call_msps.out.msp_bams,call_msps.out.msp_bams_index)                                  // create a bedgraph of 6ma calling and calculate percent 6ma coverage for each base
        pileupbedgraphtobigwig_6ma(create_pileups.out.pileups,ref_chromsizes)                                // convert the pileup bedgraph into a bigwig for downstream purposes
        pileupbedgraphtobigwig_5mC(create_pileups.out.pileups,ref_chromsizes)                                // convert the pileup bedgraph into a bigwig for downstream purposes
        pileupbedgraphtobigwig_nuc(create_pileups.out.pileups,ref_chromsizes)                                // convert the pileup bedgraph into a bigwig for downstream purposes
    }
}