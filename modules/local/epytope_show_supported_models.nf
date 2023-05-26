process EPYTOPE_SHOW_SUPPORTED_MODELS {
    label 'process_low'

    conda "bioconda::epytope=3.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.3.0--pyh7cba7a3_0' :
        'biocontainers/epytope:3.3.0--pyh7cba7a3_0' }"

    output:
    path "*.txt",        emit: txt
    path "versions.yml", emit: versions

    script:
    """
    # Extract software versions from container
    mhcflurry_version=\$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
    mhcnuggets_version=\$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")

    # Syfpeithi is not an external software, but rather a matrix on which scoring is based on -> titled version 1.0 in epytope
    syfpeithi_version=1.0

    show_supported_models.py --pred_methods syfpeithi mhcflurry mhcnuggets-class-1 mhcnuggets-class-2 --pred_method_versions \$syfpeithi_version \$mhcflurry_version \$mhcnuggets_version \$mhcnuggets_version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcflurry: \$mhcflurry_version
        mhcnuggets: \$mhcnuggets_version
        syfpeithi: \$syfpeithi_version
        epytope: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)"))
    END_VERSIONS
    """

}
