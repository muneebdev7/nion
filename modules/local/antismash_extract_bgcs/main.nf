process ANTISMASH_EXTRACT_BGCS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python_pip_biopython:9b9e6e48ba88054b' :
        'oras://community.wave.seqera.io/library/python_pip_biopython:a98a885311d28095' }"

    input:
    tuple val(meta), path(contigs), path(regions_tsv)

    output:
    tuple val(meta), path("${meta.id}_regions.fa"), emit: bgc_fasta
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('biopython'), eval('python -c "import Bio; print(Bio.__version__)"'), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample = meta.id
    """
    set -euo pipefail

    mkdir -p megahit/${sample}
    contigs_path="${contigs}"
    case "\${contigs_path}" in
        /*) ;;
        *) contigs_path="${'$'}(pwd)/\${contigs_path}" ;;
    esac
    ln -s "\${contigs_path}" "megahit/${sample}/final.contigs.fa"

    mkdir -p antismash_output/${sample}
    regions_path="${regions_tsv}"
    case "\${regions_path}" in
        /*) ;;
        *) regions_path="${'$'}(pwd)/\${regions_path}" ;;
    esac
    ln -s "\${regions_path}" "antismash_output/${sample}/${sample}_regions.tsv"

    printf "%s\n" "${sample}" > sample_list.txt

    antismash_extract_bgcs.py \\
        --sample-list sample_list.txt \\
        --megahit-dir megahit \\
        --antismash-dir antismash_output

    cp "antismash_output/${sample}/BGC_FASTAs/${sample}_regions.fa" .

    """

    stub:
    def sample = meta.id
    """
    touch ${sample}_regions.fa
    """
}
