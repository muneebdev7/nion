process METAPHLAN3_DB_DOWNLOAD {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaphlan:3.0.12--pyhb7b1952_0' :
        'biocontainers/metaphlan:3.0.12--pyhb7b1952_0' }"

    output:
    path 'metaphlan_db', emit: db

    script:
    def index_name = params.metaphlan_index ?: 'latest'
    def nproc = params.metaphlan_download_cpus ?: task.cpus
    """
    set -euo pipefail
    
    echo "[METAPHLAN3_DB_DOWNLOAD] Starting database download..."
    echo "[METAPHLAN3_DB_DOWNLOAD] Index: ${index_name}"
    echo "[METAPHLAN3_DB_DOWNLOAD] CPUs: ${nproc}"
    
    mkdir -p metaphlan_db
    
    # Download and install MetaPhlAn database
    metaphlan --install --bowtie2db metaphlan_db --index ${index_name} --nproc ${nproc} 2>&1 | tee metaphlan_install.log
    
    """

    stub:
    """
    mkdir -p metaphlan_db
    touch metaphlan_db/mpa_vJan21_CHOCOPhlAnSGB_202103.1.bt2
    touch metaphlan_db/mpa_vJan21_CHOCOPhlAnSGB_202103.rev.1.bt2
    """
}

    //
    // MODULE: Run MetaPhlAn on trimmed reads
    //
    if (params.metaphlan_db || params.metaphlan_auto_download) {
        if (params.metaphlan_db) {
            ch_metaphlan_db = channel.fromPath(params.metaphlan_db, checkIfExists: true)
        } else {
            METAPHLAN3_DB_DOWNLOAD()
            ch_metaphlan_db = METAPHLAN3_DB_DOWNLOAD.out.db
            log.warn("--metaphlan_db not provided; auto-downloading MetaPhlAn DB (index: ${params.metaphlan_index ?: 'latest'}).")
        }

        METAPHLAN3_METAPHLAN3 (
            FASTP.out.reads,
            ch_metaphlan_db
        )

        ch_versions = ch_versions.mix(METAPHLAN3_METAPHLAN3.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(METAPHLAN3_METAPHLAN3.out.profile.map { _meta, profile -> profile })
    } else {
        log.warn("MetaPhlAn is disabled: provide --metaphlan_db or set --metaphlan_auto_download true.")
    }