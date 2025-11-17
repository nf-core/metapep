process SHOW_SUPPORTED_MODELS {
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input: 
    path supported_alleles_json


    output:
    path "*.txt"       , emit: txt         
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    show_supported_models.py \\
        --json $supported_alleles_json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        supported_alleles_json: "static"
    END_VERSIONS
    """
}
