include { MERGE_FASTQ      } from '../modules/merge_fastq.nf'
include { CHOPPER          } from '../modules/chopper.nf'
include { DOWNSAMPLE_PLAN; DOWNSAMPLE } from '../modules/downsample.nf'
include { READ_STATS; READ_STATS_REPORT } from '../modules/read_stats.nf'
include { FLYE_ASSEMBLY    } from '../modules/flye.nf'
include { BANDAGE_IMAGE    } from '../modules/bandage_image.nf'
include { MEDAKA_POLISH    } from '../modules/medaka.nf'
include { PHAROKKA_DB      } from '../modules/pharokka_db.nf'
include { PHAROKKA         } from '../modules/pharokka.nf'
include { CHECKV_DB        } from '../modules/checkv_db.nf'
include { CHECKV           } from '../modules/checkv.nf'

// a size as a whole number of bases: "40kb", "0.5 Mb", "39937" -> 40000 /
// 500000 / 39937; blank -> 0 (estimate it). Stops the run at launch on
// anything else, rather than after merging and filtering every sample.
// 1 kb to 1 Mb: past the tailed-phage range (15-500 kb) at both ends, so it
// never refuses a real jumbo phage, but catches a forgotten "kb"
def sizeBp(name, value) {
    def text = (value == null || value == true) ? '' : value.toString().trim()
    if (!text) {
        return 0
    }
    def m = text =~ /^(\d+(?:\.\d+)?)\s*([kKmMgG]?)(?:[bB][pP]?)?$/
    if (!m.matches()) {
        error "--${name} '${text}' isn't a size: give a number of bases, e.g. 40kb, 170kb or 39937"
    }
    def scale = [k: 1e3, m: 1e6, g: 1e9][m.group(2).toLowerCase()] ?: 1
    def bp = (new BigDecimal(m.group(1)) * scale).toLong()
    if (bp < 1000 || bp > 1_000_000) {
        error "--${name} '${text}' is ${bp} bases; expected 1 kb to 1 Mb for a phage${bp < 1000 ? ' (did you leave off the kb?)' : ''}"
    }
    return bp
}

// Under -stub-run (CI) the database steps are skipped and PHAROKKA/CHECKV get
// tests/data/stub/{pharokka,checkv} as their `db`, which their stub blocks
// copy real T7 outputs from. Every other step has no stub block, so it still
// runs for real. Skipping the database steps (rather than stubbing them) also
// means a stub run can never leave an empty "database" in the real db_cache
def stubDb(name) {
    return channel.value(file("${projectDir}/tests/data/stub/${name}", checkIfExists: true))
}

workflow ULANA_PHAGE {

    take:
    reads_ch   // tuple(sample, fastq)

    main:
    def genome_size = sizeBp('genome_size', params.genome_size)
    if (!(params.target_coverage instanceof Integer) || params.target_coverage < 1) {
        error "--target_coverage must be a whole number of 1 or more, not '${params.target_coverage}'"
    }
    // the range LRGE's estimate is held to (modules/downsample.nf)
    if (!(params.phage_size_min instanceof Integer) || !(params.phage_size_max instanceof Integer)
            || params.phage_size_min < 1000 || params.phage_size_min >= params.phage_size_max) {
        error "--phage_size_min / --phage_size_max must be whole numbers of bases with 1000 <= min < max, not ${params.phage_size_min} / ${params.phage_size_max}"
    }

    // 1. merge multi-part fastq(.gz) files per sample
    MERGE_FASTQ(reads_ch)

    // 2. quality + length filter
    CHOPPER(MERGE_FASTQ.out.merged)

    // 2b. read length / Q-score summary, before vs. after filtering -- a
    // side branch, nothing downstream depends on it
    if (params.enable_read_stats) {
        READ_STATS(MERGE_FASTQ.out.merged.join(CHOPPER.out.filtered))
        READ_STATS_REPORT(
            READ_STATS.out.stats.map { sample, stats, hist -> stats }.collect(),
            READ_STATS.out.stats.map { sample, stats, hist -> hist }.collect()
        )
    }

    // 2c. cap the depth at --target_coverage before assembly (modules/downsample.nf);
    // a sample that doesn't reach it keeps every filtered read. Flye and Medaka
    // both get these reads
    if (params.enable_downsampling) {
        DOWNSAMPLE_PLAN(CHOPPER.out.filtered, genome_size)
        def planned = CHOPPER.out.filtered
            .join(DOWNSAMPLE_PLAN.out.plan)
            .branch { sample, reads, decision, target_bases ->
                cut: decision == 'downsample'
                keep: true
            }
        DOWNSAMPLE(planned.cut.map { sample, reads, decision, target_bases -> tuple(sample, reads, target_bases as long) })
        assembly_reads = planned.keep
            .map { sample, reads, decision, target_bases -> tuple(sample, reads) }
            .mix(DOWNSAMPLE.out.reads)
    } else {
        assembly_reads = CHOPPER.out.filtered
    }

    // 3. assemble
    FLYE_ASSEMBLY(assembly_reads)

    // 3b. render the assembly graph -- off the Flye GFA, so unaffected by Medaka
    if (params.enable_bandage) {
        BANDAGE_IMAGE(FLYE_ASSEMBLY.out.gfa)
    }

    // 4. polish -- optional; downstream steps consume whichever of these is current
    if (params.enable_medaka) {
        MEDAKA_POLISH(assembly_reads.join(FLYE_ASSEMBLY.out.assembly))
        assembly_ch = MEDAKA_POLISH.out.polished
    } else {
        assembly_ch = FLYE_ASSEMBLY.out.assembly
    }

    // 5-6 need 8.2 GB of databases; under -stub-run (CI) they're swapped for
    // the stub outputs, see stubDb() above

    // 5. annotate -- optional
    if (params.enable_pharokka) {
        if (workflow.stubRun) {
            pharokka_db = stubDb('pharokka')
        } else {
            PHAROKKA_DB()
            pharokka_db = PHAROKKA_DB.out.db
        }
        PHAROKKA(assembly_ch, pharokka_db)
    }

    // 6. completeness / contamination -- optional, independent of Pharokka
    if (params.enable_checkv) {
        if (workflow.stubRun) {
            checkv_db = stubDb('checkv')
        } else {
            CHECKV_DB()
            checkv_db = CHECKV_DB.out.db
        }
        CHECKV(assembly_ch, checkv_db)
    }

    emit:
    assembly = assembly_ch
}
