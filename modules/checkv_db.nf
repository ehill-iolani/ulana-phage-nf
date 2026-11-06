// CheckV's reference database, v1.5 (1.7 GB download, 6.4 GB unpacked once
// the DIAMOND index is built). Cached in params.db_cache by storeDir, so it's
// fetched once per cache location.
//
// Pinned to v1.5 rather than `checkv download_database`, which fetches
// whatever CURRENT_RELEASE.txt names at the time: a new release would then
// change completeness estimates between runs with no change to the pipeline.
// These are the same steps that command runs.
process CHECKV_DB {
    label 'process_low'
    container "${params.container_registry}/biocontainers/checkv:1.1.1--pyh106432d_1"
    storeDir "${params.db_cache}/checkv"

    output:
    path "checkv-db-v1.5", emit: db

    script:
    // --retry-all-errors + -C - so a dropped connection carries on from the
    // bytes already saved (see shotgun-metagenomics-nf's FETCH_DEACON_INDEX)
    """
    curl -fsSL --retry 5 --retry-all-errors --retry-delay 10 -C - -o checkv-db-v1.5.tar.gz \\
        https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
    tar -xzf checkv-db-v1.5.tar.gz
    rm checkv-db-v1.5.tar.gz

    diamond makedb --threads ${task.cpus} \\
        --in checkv-db-v1.5/genome_db/checkv_reps.faa \\
        --db checkv-db-v1.5/genome_db/checkv_reps
    """
}
