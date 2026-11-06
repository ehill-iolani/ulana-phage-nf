// Pharokka's databases (PHROGs, CARD, VFDB, the INPHARED Mash sketch): a
// 656 MB download, 1.8 GB unpacked. Cached in params.db_cache by storeDir, so
// it's fetched once per cache location, not once per run. `pharokka install`
// fetches the release that matches this Pharokka version (v1.8.0 databases
// for 1.10.1), so the directory is named after it: a version bump that needs
// new databases then fetches them instead of reusing the old ones.
process PHAROKKA_DB {
    label 'process_low'
    // --pharokka_container: see modules/pharokka.nf
    container params.pharokka_container ?: "${params.container_registry}/biocontainers/pharokka:1.10.1--pyhdfd78af_0"
    storeDir "${params.db_cache}/pharokka"

    output:
    path "pharokka_v1.8.0_databases", emit: db

    script:
    """
    pharokka install -o pharokka_v1.8.0_databases
    """
}
