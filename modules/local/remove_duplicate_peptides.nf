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
    def files =         peptides.join(' ')
    def samples =       meta.sample.join(' ')
    def conditions =    meta.conditions.collect {w -> "\"$w\""}.join(' ')
    def type =          meta.type.join(' ')
    def weights =       meta.weights.join(' ')
    def bin_basename =  meta.bin_basename.collect {w -> "\"$w\""}.join(' ')

    """
    remove_duplicate_peptides.py -i $files -s $samples -c $conditions -t $type -w $weights -b $bin_basename -o predict_peptides_${allele}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
