process CREATE_RESULTS_TABLES {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), val(taxa), path(prediction)

    output:
    tuple val(meta), path("results_*")      , emit: results
    path "versions.yml"                     , emit: versions

    script:
    """
    echo "$meta" >> results_${meta.condition}.txt
    echo "$taxa" >> results_${meta.condition}.txt
    echo $prediction >> results_${meta.condition}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
