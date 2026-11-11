// Caps each sample's read depth before assembly. Very deep runs make Flye fail:
// a ~91 kb phage at ~5,000x stopped with "No disjointigs were assembled", while
// a 1/20 subsample of the same reads assembled fine. Phage isolates are often
// sequenced that deep, so this matters more here than in ulana-nf.
//
// DOWNSAMPLE_PLAN counts the filtered bases and gets the genome size, then decides:
//   - reads above --target_coverage x genome size -> DOWNSAMPLE (Rasusa) takes
//     a random subset of about that many bases;
//   - otherwise every filtered read is used.
//
// The genome size is --genome_size when given. Otherwise LRGE estimates it from
// the reads, and the estimate is held to the tailed-phage range
// (params.phage_size_min..phage_size_max, 15-500 kb): LRGE runs high on phages
// (131 kb for 40 kb T7), and an estimate far outside the range would either
// skip downsampling or cut too deep. An estimate that fails or comes back empty
// counts as the upper bound, so a shallow sample keeps every read and a deep one
// is still capped -- unlike ulana-nf, where a failed estimate keeps every read,
// which here is the very failure this step exists to prevent.
//
// The same module as ulana-nf's apart from that limit; keep fixes in step.
// The plan is published as {sample}.downsampling.tsv either way.

process DOWNSAMPLE_PLAN {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/lrge:1.0.0--hf1b6044_0"
    publishDir(path: { "${params.outdir}/${sample}/reads" }, mode: 'copy')

    input:
    tuple val(sample), path(reads)
    val genome_size   // bp from --genome_size, or 0 to estimate it

    output:
    tuple val(sample), env('DECISION'), env('TARGET_BASES'), emit: plan
    tuple val(sample), path("${sample}.downsampling.tsv"), emit: summary

    script:
    """
    # chopper writes 4-line records, so every 4th line from the 2nd is sequence.
    # %.0f, not %d: this awk clamps %d at 2^31, so a sample over ~2.1 Gb came
    # out as -2147483648 and skipped downsampling (ulana-nf, CFT073, 2026-11-10)
    READ_BASES=\$(zcat -f ${reads} | awk 'NR % 4 == 2 { b += length(\$0) } END { printf "%.0f", b }')

    MIN_SIZE=${params.phage_size_min}
    MAX_SIZE=${params.phage_size_max}
    LRGE_ESTIMATE=''
    if [ "${genome_size}" -gt 0 ]; then
        GENOME_SIZE=${genome_size}
        SOURCE=given
    else
        # seeded so a rerun makes the same estimate -- on the same CPU count:
        # the deep barcode22 set gave 111,118 bp with 4 threads and 108,984
        # with 2 (Batch, 2026-11-06), so a machine-size change moves it ~2%
        LRGE_ESTIMATE=\$(lrge -t ${task.cpus} -s 1 ${reads} 2> lrge.log || true)
        case "\$LRGE_ESTIMATE" in
            ''|*[!0-9]*|0) LRGE_ESTIMATE=''; GENOME_SIZE=\$MAX_SIZE; SOURCE=no_estimate ;;
            *)
                if [ "\$LRGE_ESTIMATE" -lt "\$MIN_SIZE" ]; then
                    GENOME_SIZE=\$MIN_SIZE; SOURCE=estimate_raised
                elif [ "\$LRGE_ESTIMATE" -gt "\$MAX_SIZE" ]; then
                    GENOME_SIZE=\$MAX_SIZE; SOURCE=estimate_lowered
                else
                    GENOME_SIZE=\$LRGE_ESTIMATE; SOURCE=estimated
                fi
                ;;
        esac
    fi

    TARGET_BASES=\$(( ${params.target_coverage} * GENOME_SIZE ))
    if [ "\$READ_BASES" -gt "\$TARGET_BASES" ]; then
        DECISION=downsample
        ASSEMBLY_BASES=\$TARGET_BASES
    else
        DECISION=all
        ASSEMBLY_BASES=\$READ_BASES
    fi

    cov() { awk -v b="\$1" -v g="\$GENOME_SIZE" 'BEGIN { printf "%.1f", b / g }'; }
    # ulana-nf's columns, plus lrge_estimate (LRGE's own number, before the limit)
    {
        printf 'sample\\tread_bases\\tgenome_size\\tgenome_size_source\\tread_coverage\\ttarget_coverage\\treads_used\\tassembly_bases\\tassembly_coverage\\tlrge_estimate\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
            "${sample}" "\$READ_BASES" "\$GENOME_SIZE" "\$SOURCE" "\$(cov \$READ_BASES)" \\
            "${params.target_coverage}" "\$([ \$DECISION = downsample ] && echo subsample || echo all)" \\
            "\$ASSEMBLY_BASES" "\$(cov \$ASSEMBLY_BASES)" "\$LRGE_ESTIMATE"
    } > ${sample}.downsampling.tsv
    """
}

process DOWNSAMPLE {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/rasusa:5.1.0--hfa8f182_0"

    input:
    tuple val(sample), path(reads), val(target_bases)

    output:
    tuple val(sample), path("${sample}.downsampled.fastq.gz"), emit: reads

    script:
    // fixed seed: the same reads give the same subset on a rerun
    """
    rasusa reads --bases ${target_bases} --seed 1 -o ${sample}.downsampled.fastq.gz ${reads}
    """
}
