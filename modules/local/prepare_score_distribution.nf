process PREPARE_SCORE_DISTRIBUTION {
    label "process_long"
    label "process_high_memory"

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"


    input:
    path predictions
    path proteins_peptides
    path entities_proteins
    path microbiomes_entities
    path conditions
    path conditions_alleles
    path alleles

    output:
    path "prediction_scores.allele_*.tsv", emit: ch_prep_prediction_scores
    path "versions.yml"                  , emit: versions

    script:
    def chunk_size    = params.downstream_chunk_size
    def mem_log_level = params.memory_usage_log_deep ? "--mem_log_level_deep" : ""
    """
    prepare_score_distribution.py --predictions "$predictions" \\
                            --protein-peptide-occ "$proteins_peptides" \\
                            --entities-proteins-occ "$entities_proteins" \\
                            --microbiomes-entities-occ "$microbiomes_entities" \\
                            --conditions "$conditions" \\
                            --condition-allele-map "$conditions_alleles" \\
                            --alleles "$alleles" \\
                            --chunk-size $chunk_size \\
                            $mem_log_level \\
                            --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
