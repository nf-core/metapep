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
    path "*_unify_peptide_lengths.log"      , emit: log
    env unified_peptide_lengths             , emit: unified_pep_lens
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Syfpeithi is not an external software, but rather a matrix on which scoring is based on -> titled version 1.0 in epytope
    def syfpeithi_version = "1.0"
    """
    unified_peptide_lengths=\$(unify_model_lengths.py \\
                                -i $samplesheet_valid \\
                                -m $params.pred_method \\
                                -pll $params.min_pep_len \\
                                -plh $params.max_pep_len \\
                                -s "_unify_peptide_lengths")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        epytope: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)"))
        syfpeithi: $syfpeithi_version
    END_VERSIONS
    """
}
