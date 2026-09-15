#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process split_pileup_by_chr{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fire'

    input:
    tuple val(sampname), path(fdr_bed), val(ref_name), path(pileup_bed), path(pileup_bed_tbi), val(chromosome)
    output:
    tuple val(sampname), path(fdr_bed), val(ref_name), path("*.no-shuffle.${chromosome}.pileup.bed.gz"), path("*tbi"), val(chromosome)

    script:
    """
    tabix -h $pileup_bed $chromosome | bgzip > ${sampname}.FIRE_dev.no-shuffle.${chromosome}.pileup.bed.gz;
    tabix -p bed ${sampname}.FIRE_dev.no-shuffle.${chromosome}.pileup.bed.gz
    """
}


process get_chrom_sizes {
    label 'small'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple  path(ref_fasta), path(ref_fai), val(ref_name)
    output:
    tuple path(ref_fai), path("*.chrom.sizes"), val(ref_name)
    script:
    """
    cut -f 1,2 $ref_fai > ${ref_name}.chrom.sizes
    """
}

process fire_locations{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(msp_bam), val(ref_name), path(msp_bam_index)
    output:
    tuple val(sampname), path("*.FIRE_locs.bam"), val(ref_name), path("*.FIRE_locs.bam.bai")
    script:
    """
    conda run -n fiberseq-qc ft fire -t $task.cpus --min-ave-msp-size 10 --min-msp 10 \
    --skip-no-m6a $msp_bam ${sampname}.${ref_name}.FIRE_locs.bam;
    samtools index -@ $task.cpus ${sampname}.${ref_name}.FIRE_locs.bam
    """
}

process extract_unfiltered_fire_locs{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'medium'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(fire_bam), val(ref_name), path(fire_bai)
    output:
    tuple val(sampname), path("*.bed"), val(ref_name)
    script:
    """
    ft fire -t $task.cpus --extract $fire_bam | LC_ALL=C sort --parallel=$task.cpus -k1,1 -k2,2n -k3,3n -k4,4 \
    | (grep -v '^#' || true) \
    > ${sampname}.FIRE_locs.bed;
    """
}

process split_locs_by_chrom{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'medium'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(fire_locs_bed), val(ref_name), val(chrom)
    output:
    tuple val(sampname), path("*.${chrom}.FIRE_loc.bed.gz"), val(ref_name), path("*.${chrom}.FIRE_loc.bed.gz.tbi"), val(chrom)
    script:
    """
    bedextract $chrom $fire_locs_bed | bgzip -@ $task.cpus > ${sampname}.${chrom}.FIRE_loc.bed.gz;
    tabix -p bed ${sampname}.${chrom}.FIRE_loc.bed.gz
    """
}

process mosdepth{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'medium'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(fire_bam), val(ref_name), path(fire_bai), path(ref_fasta)
    output:
    tuple val(sampname), path("*coverage.bed.gz"), val(ref_name), path("*.tbi")
    script:
    """
    mosdepth -f $ref_fasta -t $task.cpus tmp ${fire_bam}
    bgzip -cd tmp.per-base.bed.gz | LC_ALL=C sort --parallel=$task.cpus -k1,1 -k2,2n -k3,3n -k4,4  \
            | bgzip -@ $task.cpus\
        > ${sampname}.coverage.bed.gz
    tabix -f -p bed ${sampname}.coverage.bed.gz
    """
}

process get_fire_cov_stats{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(cov_bed), val(ref_name), path(tbi)
    output:
    tuple val(sampname), path("*.FIRE_dev.median_coverage.txt"), val(ref_name), path("*.FIRE_dev.minimum_coverage.txt"), path("*.FIRE_dev.maximum_coverage.txt")
    script:
    """
    conda run -n FIRE python3 /opt/FIRE/scripts/cov.py --input $cov_bed --sampname $sampname
    """
}

process get_fire_locs{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(fire_bam), val(ref_name), path(fire_bai)
    output:
    tuple val(sampname), path("*.shuffle-locations.bed.gz"), val(ref_name), path("*.tbi")
    script:
    """
    ft extract $fire_bam -t $task.cpus -s --all - \
            | hck -F '#ct' -F st -F en -F fiber -F strand -F HP \
            | (grep -v "^#" || true) \
            | bgzip -@ $task.cpus \
        > ${sampname}.FIRE_dev.shuffle-locations.bed.gz
    tabix -p bed ${sampname}.FIRE_dev.shuffle-locations.bed.gz
    """
}

process coverage_filter_locs{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fire:0.20'

    input:
    // 0sampname, 1overage_bed, 2ref_name, 3coverage_bed_tbi, 4shuffle_locations_bed, 6shuffle_locations_bed_tbi, 7med_cov, 9min_cov, 10max_cov

    tuple val(sampname), path(coverage_bed), val(ref_name), path(coverage_bed_tbi), path(shuffle_locations_bed), path(shuffle_locations_bad_tbi), path(med_cov), path(min_cov), path(max_cov)
    output:
    tuple val(sampname), path("*.final_cov_filtered.bed.gz"), val(ref_name), path("*.final_cov_filtered.bed.gz.tbi")
    script:
    """
    MIN=\$(cat $min_cov)
    MAX=\$(cat $max_cov)

    bedtools intersect -header -v -f 0.2 \
            -a  $shuffle_locations_bed \
            -b <(bgzip -cd $coverage_bed | awk -v MAX="\$MAX" -v MIN="\$MIN" '\$4 <= MIN || \$4 >= MAX') \
        | bgzip -@ $task.cpus \
        > ${sampname}.final_cov_filtered.bed.gz
    tabix -p bed ${sampname}.final_cov_filtered.bed.gz
    """
}

process get_shuffled_locs{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(final_cov_filtered_locs), val(ref_name), path(final_cov_filtered_locs_tbi), path(ref_fasta), path(ref_fai)
    output:
    tuple val(sampname), path("*.fiber-locations.shuffled_filtered.bed.gz"), val(ref_name)
    script:
    """
    zcat $final_cov_filtered_locs \
    | bioawk -t '{{print \$1,\$2,\$3,\$4,\$2}}' \
    | bedtools shuffle -chrom -seed 42 \
    -i - \
    -g $ref_fai \
    |  sort -k1,1 -k2,2n -k3,3n -k4,4 \
    | bgzip -@ $task.cpus \
    > ${sampname}.fiber-locations.shuffled_filtered.bed.gz
    """
}

process get_fire_pileups_shuffled{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(fire_bam), val(ref_name), path(fire_bai),  path(shuffled_locs_bed)
    output:
    tuple val(sampname), path("*.pileup.bed.gz"), val(ref_name), path("*.pileup.bed.gz.tbi")
    script:
    """
    ft pileup $fire_bam -t $task.cpus \
            --fiber-coverage --shuffle $shuffled_locs_bed \
            --no-msp --no-nuc \
            | bgzip -@ $task.cpus \
        > ${sampname}.pileup.bed.gz
        tabix -p bed ${sampname}.pileup.bed.gz 
    """
}

process create_fdr_table{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(fire_pileup), val(ref_name), path(fire_pileup_tbi), path(median_cov), path(min_cov), path(max_cov)
    output:
    tuple val(sampname), path("*.fdr.tbl"), val(ref_name)
    script:
    """
    MIN=\$(cat $min_cov)
    MAX=\$(cat $max_cov)
    conda run -n FIRE python3 /opt/FIRE/scripts/fdr-table.v2.py \
    -v 1 $fire_pileup ${sampname}.fdr.tbl \
    --max-cov \$MAX --min-cov \$MIN 
    """
}

process get_fire_pileups_no_shuffle{
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(fire_bam), val(ref_name), path(fire_bai)
    output:
    tuple val(sampname), path("*.pileup.bed.gz"), val(ref_name), path("*.pileup.bed.gz.tbi")
    script:
    """
    ft pileup -t $task.cpus \
    --haps --fiber-coverage \
    $fire_bam \
    | bgzip -@ $task.cpus \
    > ${sampname}.FIRE_dev.no-shuffle.pileup.bed.gz;
    tabix -p bed ${sampname}.FIRE_dev.no-shuffle.pileup.bed.gz
    """
}

process make_fdr_bed {
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'large'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(fdr_table), val(ref_name), path(no_shuffle_pileup), path(pileup_tbi), val(chrom)
    output:
    tuple val(sampname), path("*.FIRE_dev.FDR.${chrom}.bed.gz"), val(ref_name), path("*.FIRE_dev.FDR.${chrom}.bed.gz.tbi"), val(chrom)

    script:
    """
    conda run -n FIRE python3 /opt/FIRE/scripts/fdr-table.py -v 1 \
    --fdr-table $fdr_table \
    $no_shuffle_pileup ${sampname}.FIRE_dev.FDR.${chrom}.bed;
    bgzip -@ 8 ${sampname}.FIRE_dev.FDR.${chrom}.bed;
    tabix -p bed ${sampname}.FIRE_dev.FDR.${chrom}.bed.gz
    """
}

process get_only_FIREs {
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'medium'
   container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(fdr_bed), val(ref_name), path(fdr_bed_tbi), val(chrom), path(fire_loc_bed), path(fire_loc_bed_tbi)
    output:
    tuple val(sampname), path("*.FIREs.${chrom}.bed.gz"), val(ref_name), path("*.FIREs.${chrom}.bed.gz.tbi"), val(chrom)
    script:
    """
    HEADER=\$(bgzip -cd ${fdr_bed} | head -n 1 || true)
    NC=\$(echo "\$HEADER" | awk '{print NF}' || true)
    FIRE_CT=\$((NC+1))
    FIRE_ST=\$((NC+2))
    FIRE_EN=\$((NC+3))
    FIRE_SIZE=\$((NC+4))
    FIRE_ID=\$((NC+5))

    OUT_HEADER=\$(printf "\$HEADER\\tpeak_chrom\\tpeak_start\\tpeak_end\\tFIRE_IDs\\tFIRE_size_mean\\tFIRE_size_ssd\\tFIRE_start_ssd\\tFIRE_end_ssd")

    zcat ${fire_loc_bed} | cut -f 1-3 | awk -v OFMT="%f" '{print \$0"\t"\$3-\$2"\t"NR}' > sam.input.tsv

    ( \\
        printf "\$OUT_HEADER\\n"; \\
        zcat ${fdr_bed} \\
            | bioawk -tc hdr '(NR==1)||(\$is_local_max=="true")' \
            | csvtk filter -tT -C '\$' -f "FDR<=0.05" \
            | csvtk filter -tT -C '\$' -f "fire_coverage>1" \
            | bioawk -tc hdr '(NR==1)||(\$fire_coverage/\$coverage>=0)' \
            | bedtools intersect -wa -wb -a - \
            -b sam.input.tsv \
            | bedtools groupby -g 1-\$NC \
            -o first,median,median,collapse,mean,sstdev,sstdev,sstdev \
            -c \$FIRE_CT,\$FIRE_ST,\$FIRE_EN,\$FIRE_ID,\$FIRE_SIZE,\$FIRE_SIZE,\$FIRE_ST,\$FIRE_EN \
    ) \\
        | hck -f 1,\$FIRE_ST,\$FIRE_EN,2-\$NC,\$FIRE_SIZE- \
        | csvtk round -tT -C '\$' -n 0 -f 2,3 \
        | bedtools sort -header -i - \
        | bgzip -@ ${task.cpus} \
        > ${sampname}.FIREs.${chrom}.bed.gz

    tabix -p bed ${sampname}.FIREs.${chrom}.bed.gz
    """
}


process merge_chromosome_beds {
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'small'
    container 'cfrasier/epi-fire:0.20'

    input:
    tuple val(sampname), path(chrom_beds), val(ref_name), path(chrom_beds_tbis)

    output:
    tuple val(sampname), path("${sampname}.FIREs.${ref_name}.bed.gz"), val(ref_name), path("${sampname}.FIREs.${ref_name}.bed.gz.tbi")

    script:
    """
    set -euo pipefail
    header=\$(zcat ${chrom_beds[0]} | head -n 1 || true)
    echo "\$header" > ${sampname}.FIREs.${ref_name}.bed
    for file in ${chrom_beds}; do
        zcat \$file | tail -n +2
    done >> ${sampname}.FIREs.${ref_name}.bed
    bgzip -@ ${task.cpus} ${sampname}.FIREs.${ref_name}.bed
    tabix -p bed ${sampname}.FIREs.${ref_name}.bed.gz
    """
}

process merge_peaks {
    publishDir "$params.outdir/6_FIRE_peaks/${sampname}", mode: 'copy'
    label 'medium'
    container 'cfrasier/epi-fire:0.21'

    input:
    tuple val(sampname), path(FIRE_bed), val(ref_name), path(med_cov), path(min_cov), path(max_cov)
    output:
    tuple val(sampname), path("*.FIRE_broad.bed.gz"), val(ref_name), path("*.FIRE_broad.bed.gz.tbi")
    script:
    """
    conda run -n FIRE python /opt/FIRE/scripts/merge_fire_peaks.py -v 1 \
                -i $FIRE_bed \
                --max-cov \$(cat ${max_cov}) \
                --min-cov \$(cat ${min_cov}) \
                --min-frac-accessible 0.0 \
                --output-file ${sampname}.FIRE_dev.FIRE_broad.bed ;
    bgzip -@ ${task.cpus} ${sampname}.FIRE_dev.FIRE_broad.bed;
    tabix -p bed ${sampname}.FIRE_dev.FIRE_broad.bed.gz
    """
}

workflow call_fire_peaks {
    take:
    fire_input_ch
    references_ch
    main:

    fire_locations(fire_input_ch)
    fire_locations.out.combine(references_ch, by: 2)
        .map { row -> tuple(row[1], row[2], row[0], row[3], row[4])}
        .set {mosdepth_input_ch}

    get_chrom_sizes(references_ch)
    get_chrom_sizes.out
        .map { row -> tuple(row[0], row[1].splitCsv(sep: '\t').collect{sublist -> sublist[0]}, row[2])}
        .transpose() // 0ref_fai, 1chr_name, 2ref_name
        .filter { ref_fai, chr, ref_name -> !(chr =~ "chrUn")}
        .filter { ref_fai, chr, ref_name -> !(chr =~ "random")}
        .filter { ref_fai, chr, ref_name -> !(chr =~ "chr[MXY]")}
        .filter { ref_fai, chr, ref_name -> !(chr =~ "chrEBV")}
        .set {chrom_sizes_ch}

    // 0 ref_name, 1 sampname, 2 fire_bam, 3 fire_bai, 4 ref_fai, 5 chrom_name
    // val(sampname), path(fire_bam), val(ref_name), path(fire_bai), val(chrom)

    extract_unfiltered_fire_locs(fire_locations.out)
    extract_unfiltered_fire_locs.out.combine(chrom_sizes_ch, by: 2)
        .map { row -> tuple(row[1], row[2], row[0], row[4])}
        // 0ref_name, 1sampname, 2fire_loc_bed, 3ref_fai, 4chrom
        // tuple val(sampname), path(fire_locs_bed), val(ref_name), val(chrom)
        .set {split_locs_by_chroms_input_ch}

    split_locs_by_chrom(split_locs_by_chroms_input_ch)
    mosdepth(mosdepth_input_ch)
    get_fire_cov_stats(mosdepth.out)
    get_fire_locs(fire_locations.out)
    mosdepth.out.combine(get_fire_locs.out, by: 0)
        .combine(get_fire_cov_stats.out, by:0)
        .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[6], row[7], row[9], row[10])}
        .set {coverage_filter_input_ch}
    // 0sampname, 1overage_bed, 2ref_name, 3coverage_bed_tbi, 4shuffle_locations_bed, 6shuffle_locations_bed, 7med_cov, 9min_cov, 10max_cov

    coverage_filter_locs(coverage_filter_input_ch)
    coverage_filter_locs.out.combine(references_ch, by: 2)
        .map { row -> tuple(row[1], row[2], row[0], row[3], row[4], row[5])}
        .set {get_shuffled_locs_input_ch}
    // 0sampname, 1coverage_bed, 2ref_name, 3coverage_bed_tbi, 4fasta, 5fasta_index

    get_shuffled_locs(get_shuffled_locs_input_ch)
    
    fire_locations.out.combine(get_shuffled_locs.out, by: 0)
        .map { row -> tuple(row[0], row[1], row[2], row[3], row[4])}
        .set {get_fire_pileups_shuffled_input_ch}
    // get_fire_pileups_shuffled_input_ch.view()
    // 0sampname, 1fire_bam, 2ref_name, 3fire_bai, 4shuffled_locs_bed

    get_fire_pileups_shuffled(get_fire_pileups_shuffled_input_ch)
    
    get_fire_pileups_shuffled.out.combine(get_fire_cov_stats.out, by: 0)
        .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[6], row[7])}
        .set {create_fdr_table_input_ch}
    // create_fdr_table_input_ch.view()
    // 0sampname, 1fire_pileup, 2ref_name, 3fire_pileup_tbi, 4median_cov, 5min_cov, 6max_cov

    create_fdr_table(create_fdr_table_input_ch)
    get_fire_pileups_no_shuffle(fire_locations.out)

    create_fdr_table.out.combine(get_fire_pileups_no_shuffle.out, by: 0)
        .map { row -> tuple(row[0], row[1], row[2], row[3], row[5])}
        .set {make_fdr_bed_input_ch}
    // 0sampname, 1fdr_table, 2ref_name, 3shuffle_pileup_bed, 4bed_tbi

    get_chrom_sizes.out
        .map { row -> tuple(row[0], row[1].splitCsv(sep: '\t').collect{sublist -> sublist[0]}, row[2])}
        .transpose() // 0ref_fai, 1chr_name, 2ref_name
        .filter { ref_fai, chr, ref_name -> !(chr =~ "chrUn")}
        .filter { ref_fai, chr, ref_name -> !(chr =~ "random")}
        .filter { ref_fai, chr, ref_name -> !(chr =~ "chr[MXY]")}
        .filter { ref_fai, chr, ref_name -> !(chr =~ "chrEBV")}
        .combine(make_fdr_bed_input_ch, by:2)
        .map { row -> tuple(row[3], row[4], row[0], row[5], row[6], row[2])}
        .set {split_pileup_by_chr_input_ch}
    // 0sampname, 1fdr_table, 2ref_name, 3shuffle_pileup_bed, 4bed_tbi, 5chromosome

    split_pileup_by_chr(split_pileup_by_chr_input_ch)

    make_fdr_bed(split_pileup_by_chr.out)
    make_fdr_bed.out.combine(split_locs_by_chrom.out, by: [0, 4])
    // 0samopname, 1chrom, 2fdr_bed, 3ref_name, 4fdr_bed_tbi, 5fire_loc_bed, 6ref_name, 7fire_loc_bed_tbi
        .map { row -> tuple(row[0], row[2], row[3], row[4], row[1], row[5], row[7])}
        .set {get_only_FIREs_input_ch}
    // val(sampname), path(fdr_bed), val(ref_name), path(fdr_bed_tbi), val(chrom), path(fire_loc_bed), path(fire_loc_bed_tbi)
    

    get_only_FIREs(get_only_FIREs_input_ch)
    get_only_FIREs.out
        .collect(flat:false)
        .flatMap()
        .groupTuple(by: 0)
        .map { row -> tuple(row[0], row[1], row[2].first(), row[3])
        //     tuple(sampname, fire_beds, ref_names, fire_tbis, chroms)
        }
        .set{merge_chrom_beds_input_ch}

    merge_chromosome_beds(merge_chrom_beds_input_ch)
     merge_chromosome_beds.out.combine(get_fire_cov_stats.out, by: 0)
        .map { row -> tuple(row[0], row[1], row[2], row[4], row[6], row[7])}
        .set {merge_peaks_input_ch}
//     // 0sampname, 1FIRE_bed, 2ref_name, 3med_cov, 4min_cov, 5max_cov

    merge_peaks(merge_peaks_input_ch)
}