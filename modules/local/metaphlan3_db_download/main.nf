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

    metaphlan --install --bowtie2db metaphlan_db --index ${index_name} --nproc ${nproc} 2>&1 | tee metaphlan_install.log
    """

    stub:
    """
    mkdir -p metaphlan_db
    touch metaphlan_db/mpa_vJan21_CHOCOPhlAnSGB_202103.1.bt2
    touch metaphlan_db/mpa_vJan21_CHOCOPhlAnSGB_202103.rev.1.bt2
    """
}

