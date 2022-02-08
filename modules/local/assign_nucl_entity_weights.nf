process ASSIGN_NUCL_ENTITY_WEIGHTS {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path conditions

    output:
    path   "conditions_entities.nucl.tsv", emit: ch_nucl_conditions_entities  // condition_id, entity_name, entity_weight
    path   "versions.yml"                , emit: versions


    script:
    """
    assign_entity_weights.py \\
        --conditions $conditions \\
        --out conditions_entities.nucl.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
