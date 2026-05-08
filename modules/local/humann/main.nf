process HUMANN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/humann:3.9--b616b2952dd13a24' :
        'oras://community.wave.seqera.io/library/humann:3.9--9df367e69764412f' }"

    input:
    tuple val(meta), path(reads), path(taxonomic_profile)
    path nucleotide_db
    path protein_db

    output:
    tuple val(meta), path("*_genefamilies.tsv") , emit: genefamilies
    tuple val(meta), path("*_pathabundance.tsv"), emit: pathabundance
    tuple val(meta), path("*_pathcoverage.tsv") , emit: pathcoverage
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single = reads instanceof Path
    def input_file = single ? "${reads}" : "${prefix}_input.fastq.gz"
    """
    # Merge paired/multiple reads into a single file if needed
    ${single ? "# Single input file — no merging required" : "cat ${reads.join(' ')} > ${input_file}"}

    humann \\
        --input ${input_file} \\
        --output humann_out \\
        --nucleotide-database ${nucleotide_db} \\
        --protein-database ${protein_db} \\
        --taxonomic-profile ${taxonomic_profile} \\
        --threads ${task.cpus} \\
        --remove-temp-output \\
        ${args}

    # Rename outputs to use the sample prefix
    mv humann_out/*_genefamilies.tsv   ${prefix}_genefamilies.tsv
    mv humann_out/*_pathabundance.tsv  ${prefix}_pathabundance.tsv
    mv humann_out/*_pathcoverage.tsv   ${prefix}_pathcoverage.tsv

    # Clean up temporary merged file to save disk space
    ${single ? "" : "rm -f ${input_file}"}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(humann --version 2>&1 | grep -oP '\\d+\\.\\d+[\\.\\d]*' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genefamilies.tsv
    touch ${prefix}_pathabundance.tsv
    touch ${prefix}_pathcoverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(humann --version 2>&1 | grep -oP '\\d+\\.\\d+[\\.\\d]*' | head -1)
    END_VERSIONS
    """
}
