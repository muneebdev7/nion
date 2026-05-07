process ANTISMASH_PROCESS_SAMPLES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python_pip_biopython:9b9e6e48ba88054b' :
        'oras://community.wave.seqera.io/library/python_pip_biopython:a98a885311d28095' }"

    input:
    tuple val(meta), path(gbk_files)

    output:
    tuple val(meta), path("${meta.id}_regions.tsv"), emit: regions
    tuple val(meta), path("${meta.id}_genes.tsv")  , emit: genes
    tuple val(meta), path("${meta.id}_summary.tsv"), emit: summary
    tuple val(meta), path("${meta.id}_merged.tsv") , emit: merged
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('biopython'), eval('python -c "import Bio; print(Bio.__version__)"'), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample = meta.id
    """
    set -euo pipefail

    mkdir -p antismash_output/${sample}
    for gbk in ${gbk_files}; do
        gbk_path="\$gbk"
        case "\$gbk_path" in
            /*) ;;
            *) gbk_path="${'$'}(pwd)/\$gbk_path" ;;
        esac
        ln -s "\$gbk_path" "antismash_output/${sample}/" || true
    done

    printf "%s\n" "${sample}" > sample_list.txt

    antismash_process_gbk.py \\
        --sample-list sample_list.txt \\
        --output-dir antismash_output
    
    cp "antismash_output/${sample}/${sample}_regions.tsv" .
    cp "antismash_output/${sample}/${sample}_genes.tsv" .
    cp "antismash_output/${sample}/${sample}_summary.tsv" .
    cp "antismash_output/${sample}/${sample}_merged.tsv" .
    """
    
    stub:
    def sample = meta.id
    """
    touch ${sample}_regions.tsv
    touch ${sample}_genes.tsv
    touch ${sample}_summary.tsv
    touch ${sample}_merged.tsv
    """
}
