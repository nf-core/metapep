process FINALIZE_CONDITION_ENTITIES {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path    entrez_conditions_entities
    path    nucl_conditions_entities
    path    microbiomes_entities
    path    entities
    path    conditions_microbiomes

    output:
    path    "conditions_entities.tsv"   , emit: ch_conditions_entities  // "condition_id", "condition_name", "entity_id", "entity_weight"
    path    "versions.yml"              , emit: versions

    script:

    """
    finalize_condition_entities.py \\
        -ece $entrez_conditions_entities \\
        -nce $nucl_conditions_entities \\
        -men $microbiomes_entities \\
        -ent "$entities" \\
        -cond "$conditions_microbiomes" \\
        -o conditions_entities.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
