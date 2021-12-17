process COLLECT_STATS {
    label 'process_long'
    label 'process_high_memory'
    label 'cache_lenient'

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
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
    collect_stats.py --peptides "$peptides" \\
                    --protein-peptide-occ "$proteins_peptides" \\
                    --entities-proteins-occ "$entities_proteins" \\
                    --microbiomes-entities-no-weights-occ "$microbiomes_entities" \\
                    --conditions "$conditions" \\
                    --outfile stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """

}
