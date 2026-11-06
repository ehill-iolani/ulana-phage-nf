// Phage annotation with Pharokka: PHANOTATE gene calls, PHROG functions and
// categories, CARD/VFDB flags, tRNAs/tmRNAs/CRISPRs, and the closest phage in
// INPHARED by Mash distance.
//
// Meta mode (-m) when the assembly has more than one contig: an isolate's
// assembly can hold more than the phage (a host plasmid or fragment), and meta
// mode gives each contig its own Mash hit and stats instead of one hit for the
// file. Pharokka refuses -m on a single contig, so it's only added then. Meta
// mode switches the default gene caller to prodigal-gv and drops the HMM
// search, so -g phanotate (always) and --meta_hmm put back the isolate
// defaults. No --dnaapler: contigs keep Flye's start (decided 2026-11-06).
//
// A sample's novel phage has no INPHARED hit within --mash_distance 0.2; the
// mash file then holds the literal "no_inphared_mash_hit". That's a normal
// result, not a failure.
process PHAROKKA {
    tag "$sample"
    label 'process_medium'
    // --pharokka_container swaps the image: the stock one segfaults under
    // Docker's amd64 emulation on Apple silicon (its polars runtime needs
    // AVX2). docker/pharokka-compat/Dockerfile builds a local fix; Batch VMs
    // are real x86 and use the stock image
    container params.pharokka_container ?: "${params.container_registry}/biocontainers/pharokka:1.10.1--pyhdfd78af_0"
    publishDir(path: { "${params.outdir}/${sample}/annotation/pharokka" }, mode: 'copy',
               saveAs: { name -> name.endsWith('.version') ? null : name })

    input:
    tuple val(sample), path(assembly)
    path db

    output:
    tuple val(sample), path("pharokka.gff"), emit: gff
    tuple val(sample), path("pharokka_cds_final_merged_output.tsv"), emit: cds
    tuple val(sample), path("pharokka_top_hits_mash_inphared.tsv"), emit: closest
    tuple val(sample), path("pharokka.version"), emit: version
    tuple val(sample), path("{pharokka*,top_hits_*,phanotate.f*,terL.f*}"), emit: files

    script:
    // a fixed locus tag (not Pharokka's random one), so a rerun names genes
    // the same way; letters and digits only
    def locus = sample.replaceAll(/[^A-Za-z0-9]/, '').toUpperCase() ?: 'PHAGE'
    """
    META=''
    if [ "\$(grep -c '^>' ${assembly})" -gt 1 ]; then
        META='-m --meta_hmm'
    fi
    pharokka run -i ${assembly} -o pharokka_out -d ${db} -t ${task.cpus} -f \\
        -p pharokka -l ${locus} -g phanotate \$META

    # the results sit beside Pharokka's working dirs; bring up just the files,
    # then rename what's left so the output glob doesn't publish its logs/
    for f in pharokka_out/*; do
        if [ -f "\$f" ]; then mv "\$f" .; fi
    done
    mv pharokka_out work_pharokka
    # "pharokka_1793996051.0661829.log" -> a name the results page can find
    mv pharokka_*.log pharokka.log

    # "pharokka v1.10.1" -> "1.10.1"
    pharokka -V 2>&1 | awk '{ sub(/^v/, "", \$NF); print \$NF }' > pharokka.version
    """

    stub:
    // -stub-run (CI): the workflow passes tests/data/stub/pharokka as `db`
    // instead of fetching the databases, and this copies the real T7 outputs
    // in from it, whatever the assembly
    """
    cp -L ${db}/* .
    echo 1.10.1 > pharokka.version
    """
}
