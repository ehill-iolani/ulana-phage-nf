// Assembles each sample with Flye.
//
// Second try on very deep reads: DOWNSAMPLE_PLAN caps the depth before this
// step, but only as well as it knows the genome size. When a sample carries
// host DNA (a whole-plaque WGS prep), LRGE sizes the host, not the phage, so
// nothing is cut: WGS_Lewando (PRJNA1338505, 2026-11-10) was estimated at
// 3.7 Mb, held to 500 kb, kept all 46 Mb of reads, and Flye stopped with "No
// disjointigs were assembled" at an overlap-based coverage of 583x. Flye
// measures that coverage itself before it gives up, so when it fails that way
// at more than --target_coverage, this takes a seeded random subset of
// target/measured of the reads (150/583 = 26% there) and runs Flye once more.
// Any other failure, or a second failure, still fails the task. Only with
// --enable_downsampling, since turning that off asks for every read.
//
// Medaka still polishes with every read this step was given. flye_reads.tsv
// records which reads were assembled.
//
// One copy per genome: Flye sometimes reports a phage as one contig holding
// two back-to-back copies (8NZ0, EGVZ, ZTYS in PRJNA1502461, 2026-11-11:
// ~85 kb for ~43 kb genomes). bin/collapse_copies.py cuts such a contig to one
// copy before Medaka, so polishing, CheckV and Pharokka see the real genome
// (polishing the ZTYS dimer as it was left 368 differences; one copy, 5 plus
// one ~190 bp stretch). copies.tsv lists what was cut; assembly_info.txt and
// the graph stay as Flye made them. Contigs over 2 x --phage_size_max (host
// DNA) aren't checked.

process FLYE_ASSEMBLY {
    tag "$sample"
    label 'process_high'
    container "${params.container_registry}/biocontainers/flye:2.9.6--py312h734f728_1"
    // the .version is for the GenBank assembly comment, not a result
    publishDir(path: { "${params.outdir}/${sample}/assembly/flye" }, mode: 'copy', saveAs: { it.endsWith('.version') ? null : it })

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("assembly.fasta"), emit: assembly
    tuple val(sample), path("assembly_info.txt"), emit: info
    tuple val(sample), path("assembly_graph.gfa"), emit: gfa
    tuple val(sample), path("flye_reads.tsv"), emit: reads_used
    tuple val(sample), path("copies.tsv"), emit: copies
    tuple val(sample), path("flye.version"), emit: version

    script:
    """
    run_flye() {
        rm -rf flye_run
        flye ${params.flye_mode} "\$1" --out-dir flye_run --threads ${task.cpus}
    }
    # Flye's own depth estimate from read overlaps ("Overlap-based coverage: 583")
    overlap_cov() {
        grep -o 'Overlap-based coverage: [0-9]*' flye_run/flye.log | tail -n 1 | grep -o '[0-9]*\$' || true
    }

    READS_USED=all
    FRACTION=1
    FIRST_COV=''
    if ! run_flye ${reads}; then
        COV=\$(overlap_cov)
        if [ "${params.enable_downsampling}" = true ] \\
            && grep -q 'No disjointigs were assembled' flye_run/flye.log \\
            && [ -n "\$COV" ] && [ "\$COV" -gt ${params.target_coverage} ]; then
            FRACTION=\$(awk -v t=${params.target_coverage} -v c="\$COV" 'BEGIN { printf "%.4f", t / c }')
            echo "Flye assembled no disjointigs at \${COV}x; retrying on a seeded \$FRACTION subset of the reads" >&2
            FIRST_COV=\$COV
            # 4-line fastq records (chopper/Rasusa output); seed 1 so a rerun
            # picks the same reads
            python3 - ${reads} "\$FRACTION" <<'PY'
import gzip, random, sys
src, frac = sys.argv[1], float(sys.argv[2])
rng = random.Random(1)
opener = gzip.open if src.endswith('.gz') else open
with opener(src, 'rb') as fin, gzip.open('flye_subset.fastq.gz', 'wb', compresslevel=1) as fout:
    while True:
        rec = [fin.readline() for _ in range(4)]
        if not rec[0]:
            break
        if rng.random() < frac:
            fout.writelines(rec)
PY
            READS_USED=subset
            run_flye flye_subset.fastq.gz
        else
            exit 1
        fi
    fi
    mv flye_run/assembly_info.txt flye_run/assembly_graph.gfa .

    # a contig that is two copies of the genome -> one copy (see above)
    collapse_copies.py --sample "${sample}" --in flye_run/assembly.fasta \\
        --out assembly.fasta --table copies.tsv --max-length ${(params.phage_size_max as long) * 2}

    # first_try_coverage = Flye's overlap coverage on the reads it was given,
    # assembled_coverage = the same on the reads it assembled
    {
        printf 'sample\\treads_used\\tsubset_fraction\\tfirst_try_coverage\\tassembled_coverage\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${sample}" "\$READS_USED" "\$FRACTION" \\
            "\${FIRST_COV:-\$(overlap_cov)}" "\$(overlap_cov)"
    } > flye_reads.tsv

    # "2.9.6-b1802" -> "2.9.6"
    flye --version | cut -d- -f1 > flye.version
    """
}
