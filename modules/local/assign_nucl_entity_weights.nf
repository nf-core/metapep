process ASSIGN_NUCL_ENTITY_WEIGHTS {
    tag "$microbiome_ids"
    label 'process_low'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    val  microbiome_ids
    path weights_files

    output:
    path "microbiomes_entities.nucl.tsv", emit: ch_nucl_microbiomes_entities  // entity_name, microbiome_id, entity_weight
    path "versions.yml"                 , emit: versions


    script:
    microbiome_ids = microbiome_ids.join(' ')
    """
    assign_entity_weights.py \\
        --microbiome-ids $microbiome_ids \\
        --weights-files $weights_files \\
        --out microbiomes_entities.nucl.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
