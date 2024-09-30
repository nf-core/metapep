process FINALIZE_MICROBIOME_ENTITIES {
    label 'process_low'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path(entrez_microbiomes_entities)
    path(nucl_microbiomes_entities)
    path(microbiomes_entities_noweights)
    path(entities)

    output:
    path    "microbiomes_entities.tsv"  , emit: ch_microbiomes_entities  // entity_id, microbiome_id, entity_weight
    path    "versions.yml"              , emit: versions

    script:

    """
    finalize_microbiome_entities.py \\
        -eme $entrez_microbiomes_entities \\
        -nme $nucl_microbiomes_entities \\
        -menw $microbiomes_entities_noweights \\
        -ent "$entities" \\
        -o microbiomes_entities.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
