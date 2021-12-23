process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path samplesheet

    output:
    path "microbiomes.csv"          , emit: microbiomes                  // microbiome_id, microbiome_path, microbiome_type, weights_path
    path "conditions.csv"           , emit: conditions                   // condition_id, condition_name, microbiome_id
    path "alleles.csv"              , emit: alleles                      // allele_id, allele_name
    path "conditions_alleles.csv"   , emit: conditions_alleles           // condition_id, allele_id
    path "versions.yml"             , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/metapep/bin/
    """
    check_samplesheet.py \\
        -i $samplesheet \
        -m microbiomes.csv \
        -c conditions.csv \
        -a alleles.csv \
        -ca conditions_alleles.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
