process PREPARE_ENTITY_BINDING_RATIOS {
    label "process_long"
    label "process_high_memory"

    conda (params.enable_conda ? "conda-forge::pandas=1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    path predictions
    path proteins_peptides
    path entities_proteins
    path microbiomes_entities
    path conditions
    path conditions_alleles
    path alleles

    output:
    path "entity_binding_ratios.allele_*.tsv"   , emit: ch_prep_entity_binding_ratios
    path "versions.yml"                         , emit: versions

    script:
    def proc_chunk_size       = params.proc_chunk_size  // TODO add independent parameter oand/r document
    """
    prepare_entity_binding_ratios.py --predictions "$predictions" \\
                            --protein-peptide-occ "$proteins_peptides" \\
                            --entities-proteins-occ "$entities_proteins" \\
                            --microbiomes-entities-occ "$microbiomes_entities" \\
                            --conditions "$conditions" \\
                            --condition-allele-map "$conditions_alleles" \\
                            --alleles "$alleles" \\
                            --method ${params.pred_method} \\
                            --chunk-size $proc_chunk_size \\
                            --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
