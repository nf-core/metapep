process PREPARE_SCORE_DISTRIBUTION {
    label "process_long"
    label "process_high_memory"

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path predictions
    path proteins_peptides
    path entities_proteins
    path conditions_entities
    path conditions
    path conditions_alleles
    path alleles

    output:
    path "prediction_scores.allele_*.tsv", emit: ch_prep_prediction_scores
    path "versions.yml"                  , emit: versions

    script:
    """
    prepare_score_distribution.py --predictions "$predictions" \\
                            --protein-peptide-occ "$proteins_peptides" \\
                            --entities-proteins-occ "$entities_proteins" \\
                            --conditions-entities-occ "$conditions_entities" \\
                            --conditions "$conditions" \\
                            --condition-allele-map "$conditions_alleles" \\
                            --alleles "$alleles" \\
                            --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
