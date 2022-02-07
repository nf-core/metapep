process SORT_PEPTIDES {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_sorted.tsv"), emit: output

    script:
    def name = input.getBaseName()

    """
    (head -n 1 $input && tail -n +2 $input | sort) > "${name}_sorted.tsv"
    """
}
