process PREPARE_PLOTS {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(predictions)

    output:
    tuple val(meta), path("prediction_scores.*.tsv"),         emit: ch_prep_prediction_scores
    tuple val(meta), path("entity_binding_ratios.*.tsv"),     emit: ch_prep_entity_binding_ratios
    path "versions.yml",                    emit: versions


    script:

    """
    prepare_plots.py -i $predictions -a $meta.alleles

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
