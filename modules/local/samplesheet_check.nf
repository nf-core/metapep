process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path samplesheet

    output:
    path "microbiomes.tsv"          , emit: microbiomes                  // microbiome_id, microbiome_path, microbiome_type, weights_path
    path "conditions.tsv"           , emit: conditions                   // condition_id, condition_name, microbiome_id
    path "alleles.tsv"              , emit: alleles                      // allele_id, allele_name
    path "conditions_alleles.tsv"   , emit: conditions_alleles           // condition_id, allele_id
    path "samplesheet.valid.csv"    , emit: samplesheet_valid
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/metapep/bin/
    """
    check_samplesheet.py \\
        -i $samplesheet \\
        -m microbiomes.tsv \\
        -c conditions.tsv \\
        -a alleles.tsv \\
        -ca conditions_alleles.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
