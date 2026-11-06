// Completeness and host contamination per contig with CheckV, replacing
// ulana-nf's CheckM (which has no phage marker sets: a complete phage scores
// 0% there). Every contig is kept and labelled; quality_summary.tsv says which
// look viral and which look like host.
//
// Flye calls a terminally redundant linear phage like T7 "circular" (its ends
// overlap); CheckV's complete_genomes.tsv only lists contigs with a direct or
// inverted terminal repeat it can see in the sequence, so expect it to be
// empty for most Flye assemblies.
process CHECKV {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/checkv:1.1.1--pyh106432d_1"
    publishDir(path: { "${params.outdir}/${sample}/qc/checkv" }, mode: 'copy')

    input:
    tuple val(sample), path(assembly)
    path db

    output:
    tuple val(sample), path("quality_summary.tsv"), emit: summary
    tuple val(sample), path("{completeness,contamination,complete_genomes}.tsv"), emit: details

    script:
    """
    checkv end_to_end ${assembly} checkv_out -d ${db} -t ${task.cpus} --remove_tmp
    mv checkv_out/{quality_summary,completeness,contamination,complete_genomes}.tsv .
    """

    stub:
    // -stub-run (CI): the workflow passes tests/data/stub/checkv as `db`
    // instead of fetching the database, and this copies the real T7 outputs
    // in from it, whatever the assembly
    """
    cp -L ${db}/* .
    """
}
