process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "bioconda::epytope=3.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.3.0--pyh7cba7a3_0' :
        'quay.io/biocontainers/epytope:3.3.0--pyh7cba7a3_0' }"

    input:
    path samplesheet

    output:
    path "microbiomes.tsv"          , emit: microbiomes                  // microbiome_id, microbiome_path, microbiome_type, weights_path, microbiome_bare_id
    path "conditions.tsv"           , emit: conditions                   // condition_id, condition_name, microbiome_id
    path "alleles.tsv"              , emit: alleles                      // allele_id, allele_name
    path "conditions_alleles.tsv"   , emit: conditions_alleles           // condition_id, allele_id
    path "samplesheet.valid.csv"    , emit: samplesheet_valid
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/metapep/bin/
    def pred_method           = params.pred_method
    switch (pred_method) {
        case "syfpeithi":
            pred_method_version = "1.0";
            break;
        case "mhcflurry":
            pred_method_version = "1.4.3";
            break;
        case "mhcnuggets-class-1":
            pred_method_version = "2.3.2";
            break;
        case "mhcnuggets-class-2":
            pred_method_version = "2.3.2";
            break;
        }
    """
    check_samplesheet.py \\
        -i $samplesheet \\
        -m microbiomes.tsv \\
        -c conditions.tsv \\
        -a alleles.tsv \\
        -ca conditions_alleles.tsv \\
        -pm ${pred_method} \\
        -pmv ${pred_method_version}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
