process SPLIT_PRED_TASKS {
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
    path(conditions_alleles  )
    path(alleles             )
    // The tables are joined to map peptide -> protein -> microbiome -> condition -> allele
    // and thus to enumerate, which (peptide, allele) combinations have to be predicted.

    output:
    path "peptides_*.txt",  emit:   ch_epitope_prediction_chunks
    path "versions.yml",    emit:   versions

    script:
    def pred_chunk_size       = params.pred_chunk_size
    def subsampling = params.sample_n > 0 ? "--sample_n ${params.sample_n}" : ""
    """
    gen_prediction_chunks.py --peptides "$peptides" \\
                            --protein-peptide-occ "$proteins_peptides" \\
                            --entities-proteins-occ "$entities_proteins" \\
                            --microbiomes-entities-occ "$microbiomes_entities" \\
                            --conditions "$conditions" \\
                            --condition-allele-map "$conditions_alleles" \\
                            --max-chunk-size $pred_chunk_size \\
                            $subsampling \\
                            --alleles "$alleles" \\
                            --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """

}
