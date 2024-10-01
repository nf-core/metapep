process CHECK_SAMPLESHEET_CREATE_TABLES {
    tag "$samplesheet"
    label 'process_single'

    conda "bioconda::epytope=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.3.1--pyh7cba7a3_0' :
        'biocontainers/epytope:3.3.1--pyh7cba7a3_0' }"

    input:
    path samplesheet

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
    # Extract software versions from container
    mhcflurry_version=\$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
    mhcnuggets_version=\$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")

    # Syfpeithi is not an external software, but rather a matrix on which scoring is based on -> titled version 1.0 in epytope
    syfpeithi_version=1.0

    # Assign version based on method
    case $params.pred_method in

        "syfpeithi")
        pred_method_version=\$syfpeithi_version
        ;;

        "mhcflurry")
        pred_method_version=\$mhcflurry_version
        ;;

        "mhcnuggets-class-1" | "mhcnuggets-class-2")
        pred_method_version=\$mhcnuggets_version
        ;;
    esac


    check_samplesheet_create_tables.py \\
        -i $samplesheet \\
        -m microbiomes.tsv \\
        -c conditions.tsv \\
        -a alleles.tsv \\
        -ca conditions_alleles.tsv \\
        -pm $params.pred_method \\
        -pmv \$pred_method_version \\
        -pl $params.min_pep_len $params.max_pep_len

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        epytope: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)"))
        mhcflurry: \$mhcflurry_version
        mhcnuggets: \$mhcnuggets_version
        syfpeithi: \$syfpeithi_version
    END_VERSIONS
    """
}
