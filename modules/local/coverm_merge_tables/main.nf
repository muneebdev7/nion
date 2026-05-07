process COVERM_MERGE_TABLES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python_pip_biopython:9b9e6e48ba88054b' :
        'oras://community.wave.seqera.io/library/python_pip_biopython:a98a885311d28095' }"

    input:
    tuple val(meta), path(coverage_tables)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: merged
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "coverm_merged"
    def input_list = coverage_tables.sort { path_item -> path_item.toString() }.join(" ")
    """
    coverm_sanitize_and_merge.py \\
        ${args} \\
        --output ${prefix}.tsv \\
        ${input_list}

    """

    stub:
    prefix = task.ext.prefix ?: "coverm_merged"
    """
    touch ${prefix}.tsv

    """
}
