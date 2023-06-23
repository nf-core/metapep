process PREPARE_SCORE_DISTRIBUTION {
    label "process_long"
    label "process_high_memory"

    // TODO generate extra biocontainer with only specifying pandas version (currently mulled container taken from "bedtools=2.23.0,pandas=1.5.2")
    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d19e2715c83e4582e3f1fb0a2e473abde8ca636e:fc171b36fc2e2a38a259a1c82a139b59d94c968b-0' :
        'biocontainers/mulled-v2-d19e2715c83e4582e3f1fb0a2e473abde8ca636e:fc171b36fc2e2a38a259a1c82a139b59d94c968b-0' }"

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
    def chunk_size            = params.chunk_size * params.chunk_size_scaling
    def mem_log_level         = params.memory_usage_log_deep ? "--mem_log_level_deep" : ""
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
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
