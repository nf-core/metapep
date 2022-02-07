process REMOVE_DUPLICATE_PEPTIDES  {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(allele), val(meta), val(peptides)

    output:
    tuple val(allele), path("predict_peptides_*.tsv"), emit: output
    path "versions.yml", emit: versions

    script:
    def files = peptides.join(' ')
    def entities = meta.sample.join(' ')
    """
    remove_duplicate_peptides.py -i $files -e $entities -o predict_peptides_${allele}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
