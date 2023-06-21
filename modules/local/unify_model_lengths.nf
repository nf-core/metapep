process UNIFY_MODEL_LENGTHS {
    tag "$samplesheet_valid"
    label 'process_single'

    conda "bioconda::epytope=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.3.1--pyh7cba7a3_0' :
        'biocontainers/epytope:3.3.1--pyh7cba7a3_0' }"

    input:
    path samplesheet_valid

    output:
    path "unified_allele_models.tsv"    , emit: allele_models  // allele_name, peptide_length, allele_model
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Syfpeithi is not an external software, but rather a matrix on which scoring is based on -> titled version 1.0 in epytope
    def syfpeithi_version = "1.0"
    """
    unify_model_lengths.py \\
        -i $samplesheet_valid \\
        -m $params.pred_method \\
        -pll $params.min_pep_len \\
        -plh $params.max_pep_len \\
        -o "unified_allele_models.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        epytope: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)"))
        mhcflurry: \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
        mhcnuggets: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")
        syfpeithi: $syfpeithi_version
    END_VERSIONS
    """
}