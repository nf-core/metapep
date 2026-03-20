process CHECK_SAMPLESHEET_CREATE_TABLES {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::mhcgnomes:1.8.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.6--pyh7cba7a3_0' :
        'biocontainers/mhcgnomes:1.8.6--pyh7cba7a3_0' }"

    input:
    path samplesheet
   // path supported_alleles_json

    output:
    path "microbiomes.tsv"       , emit: microbiomes                  // microbiome_id, microbiome_path, microbiome_type, weights_path, microbiome_bare_id
    path "conditions.tsv"        , emit: conditions                   // condition_id, condition_name, microbiome_id
    path "alleles.tsv"           , emit: alleles                      // allele_id, allele_name
    path "conditions_alleles.tsv", emit: conditions_alleles           // condition_id, allele_id
    path "samplesheet.valid.csv" , emit: samplesheet_valid
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_samplesheet_create_tables.py \\
        -i $samplesheet \\
        -m microbiomes.tsv \\
        -c conditions.tsv \\
        -a alleles.tsv \\
        -ca conditions_alleles.tsv \\
        -pm $params.pred_method \\
        -pl $params.min_pep_len $params.max_pep_len \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        mhcgnomes: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcgnomes').version)")
    END_VERSIONS
    """
}
