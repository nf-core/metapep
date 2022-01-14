process COLLECT_STATS {
    label 'process_long'
    label 'process_high_memory'

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"
    
    input:
    path(peptides            )
    path(proteins_peptides   )
    path(entities_proteins   )
    path(microbiomes_entities)
    path(conditions          )

    output:
    path "stats.txt",       emit:   ch_stats
    path "versions.yml",    emit:   versions

    script:
    """
    collect_stats.py --peptides "$peptides" \
                     --protein-peptide-occ "$proteins_peptides" \
                     --entities-proteins-occ "$entities_proteins" \
                     --microbiomes-entities-occ "$microbiomes_entities" \
                     --conditions "$conditions" \
                     --outfile stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}