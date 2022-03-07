process REMOVE_DUPLICATE_PEPTIDES  {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(allele), val(meta), path(peptides)
    path(weights_table)

    output:
    tuple val(allele), path("predict_peptides_*.tsv.gz"), emit: output
    path "versions.yml", emit: versions

    script:
    def samples =       meta.sample.join(' ')
    def conditions =    meta.conditions.collect {w -> "\"$w\""}.join(' ')
    def type =          meta.type.join(' ')
    def weights_ids =   meta.weights.collect {w -> "\"$w\"" ?: "null"}.join(' ')
    def bin_basename =  meta.bin_basename.collect {w -> "\"$w\""}.join(' ')

    """
    remove_duplicate_peptides.py -i $peptides -s $samples -c $conditions -t $type -wi $weights_ids -w $weights_table -b $bin_basename -f -o predict_peptides_${allele}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
