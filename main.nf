#!/usr/bin/env nextflow
/*
 * ULANA-PHAGE-NF -- phage genome assembly and annotation from Oxford Nanopore
 * long reads. Built from ulana-nf's skeleton: the same read handling, Flye
 * assembly and Medaka polishing (those modules are copies; keep fixes in
 * step), with ulana-nf's bacterial steps (Prokka, CheckM, identifier genes,
 * AMRFinder, GenBank) swapped for Pharokka annotation and CheckV QC.
 */

nextflow.enable.dsl = 2

include { ULANA_PHAGE } from './workflows/ulana_phage.nf'
include { MERGE_FASTQ } from './modules/merge_fastq.nf'
include { CHOPPER     } from './modules/chopper.nf'
include { READ_STATS; READ_STATS_REPORT } from './modules/read_stats.nf'

// ---- top-level params (override via -params-file or --flag) ----
// tool params and the enable_* toggles all live in nextflow.config, not here
// -- anything read directly inside an included module or workflow script
// must be set there to be reliably visible by the time that script runs
params.input  = null   // path to samplesheet.csv (sample,fastq); fastq globs resolve relative to the launch dir, not the CSV's location
params.outdir = "results"
params.help   = false

def helpMessage() {
    log.info """
    ULANA-PHAGE-NF -- phage genome assembly and annotation (ONT)
    ------------------------------------------------------------------
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results

    Required:
      --input       CSV: sample,fastq (fastq may be a single file or a glob)

    Optional steps:
      --enable_read_stats    read length/Q-score summary before vs. after filtering (default ${params.enable_read_stats})
      --enable_downsampling  subsample very deep samples to --target_coverage before assembly (default ${params.enable_downsampling})
      --enable_medaka        polish the Flye assembly with Medaka (default ${params.enable_medaka})
      --enable_bandage       render a PNG of the Flye assembly graph (default ${params.enable_bandage})
      --enable_pharokka      annotate with Pharokka: PHROG functions, closest INPHARED phage, CARD/VFDB (default ${params.enable_pharokka})
      --enable_checkv        completeness and host contamination with CheckV (default ${params.enable_checkv})

    Key tool params:
      --chopper_q / --chopper_minlength   chopper quality/length filtering (default ${params.chopper_q} / ${params.chopper_minlength})
      --target_coverage                   depth to subsample to; a sample below it keeps every read (default ${params.target_coverage})
      --genome_size                       genome size for --target_coverage, e.g. 40kb (default: LRGE estimates it per
                                          sample, held to ${params.phage_size_min}-${params.phage_size_max} bp)
      --flye_mode                         Flye read-type flag (default ${params.flye_mode})
      --medaka_model                      Medaka model name (default ${params.medaka_model})
      --pharokka_container                Pharokka image to use instead of the biocontainer (Apple-silicon Macs, see README)

    See nextflow_schema.json for the full parameter list.
    """.stripIndent()
}

// Two rows with the same sample name write the same output files, and
// Nextflow only notices at the first step that gathers every sample. Read the
// sheet up front so the run fails before anything is submitted. Trimmed and
// case-insensitive, like the platform's submit form.
def checkUniqueSamples(path) {
    def names = file(path, checkIfExists: true).splitCsv(header: true).collect { (it.sample ?: '').trim() }
    def repeated = names.groupBy { it.toLowerCase() }.findAll { key, group -> key && group.size() > 1 }
    if (repeated) {
        error "Samplesheet: ${repeated.values()*.first().collect { "'${it}'" }.join(', ')} used by more than one row. Each sample's outputs are named after it, so every row needs its own sample name (list one sample's several files in a single row with a glob)"
    }
}

def samplesheetToChannel(path) {
    checkUniqueSamples(path)
    Channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // checkIfExists doesn't catch a glob that matches zero files, so check explicitly --
            // otherwise this silently produces an empty file list downstream instead of failing
            def fq = file(row.fastq, checkIfExists: true)
            def files = fq instanceof List ? fq : [fq]
            if (files.isEmpty()) {
                error "Samplesheet row '${row.sample}': no files matched '${row.fastq}' (glob paths resolve relative to the launch directory, not the samplesheet's location)"
            }
            tuple(row.sample, fq)
        }
}

workflow {
    if (params.help || !params.input) {
        helpMessage()
        exit 0
    }

    ULANA_PHAGE(samplesheetToChannel(params.input))
}

// dev/debug entry point -- run just the merge step in isolation, e.g.:
//   nextflow run main.nf -entry MERGE_ONLY --input samplesheet.csv -profile docker
workflow MERGE_ONLY {
    if (!params.input) {
        log.error "MERGE_ONLY requires --input samplesheet.csv"
        exit 1
    }
    MERGE_FASTQ(samplesheetToChannel(params.input))
}

// Read QC only: merge -> chopper -> before/after length & Q-score summary, no
// assembly. Handy for picking --chopper_q/--chopper_minlength first:
//   nextflow run main.nf -entry READ_QC_ONLY --input samplesheet.csv --chopper_q 12 -profile docker
workflow READ_QC_ONLY {
    if (!params.input) {
        log.error "READ_QC_ONLY requires --input samplesheet.csv"
        exit 1
    }
    MERGE_FASTQ(samplesheetToChannel(params.input))
    CHOPPER(MERGE_FASTQ.out.merged)
    READ_STATS(MERGE_FASTQ.out.merged.join(CHOPPER.out.filtered))
    READ_STATS_REPORT(
        READ_STATS.out.stats.map { sample, stats, hist -> stats }.collect(),
        READ_STATS.out.stats.map { sample, stats, hist -> hist }.collect()
    )
}
