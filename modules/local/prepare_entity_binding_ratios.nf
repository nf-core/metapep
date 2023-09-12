process PREPARE_ENTITY_BINDING_RATIOS {
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
    path "entity_binding_ratios.allele_*.tsv", emit: ch_prep_entity_binding_ratios
    path "versions.yml"                      , emit: versions

    script:
    def chunk_size                = params.downstream_chunk_size
    def syfpeithi_score_threshold = params.syfpeithi_score_threshold
    def mhcf_mhcn_score_threshold = params.mhcflurry_mhcnuggets_score_threshold
    def mem_log_level             = params.memory_usage_log_deep ? "--mem_log_level_deep" : ""
    """
    prepare_entity_binding_ratios.py --predictions "$predictions" \\
                            --protein-peptide-occ "$proteins_peptides" \\
                            --entities-proteins-occ "$entities_proteins" \\
                            --microbiomes-entities-occ "$microbiomes_entities" \\
                            --conditions "$conditions" \\
                            --condition-allele-map "$conditions_alleles" \\
                            --alleles "$alleles" \\
                            --method ${params.pred_method} \\
                            --chunk-size $chunk_size \\
                            --syfpeithi_score_threshold $syfpeithi_score_threshold \\
                            --mhcf_mhcn_score_threshold $mhcf_mhcn_score_threshold \\
                            $mem_log_level \\
                            --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
