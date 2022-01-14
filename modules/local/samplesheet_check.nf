process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path samplesheet

    output:
    path "microbiomes.tsv"          , emit: microbiomes                  // microbiome_id, microbiome_path, microbiome_type, weights_path
    path "conditions.tsv"           , emit: conditions                   // condition_id, condition_name, microbiome_id
    path "alleles.tsv"              , emit: alleles                      // allele_id, allele_name
    path "conditions_alleles.tsv"   , emit: conditions_alleles           // condition_id, allele_id
    path "versions.yml"             , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/metapep/bin/
    """
    check_samplesheet.py \\
        -i $samplesheet \
        -m microbiomes.tsv \
        -c conditions.tsv \
        -a alleles.tsv \
        -ca conditions_alleles.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
