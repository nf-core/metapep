process MERGE_CHUNKS_BUFFER {

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path    predictions
    path    PEPTIDE_MAP
    path    ALLELE_MAP

    output:
    path "predictions.buffer_*.tsv", emit: ch_predictions_merged_buffer
    path "versions.yml"            , emit: versions

    script:
    def chunk_size = params.prediction_chunk_size * params.pred_chunk_size_scaling
    """
     [[ ${predictions[0]} =~ peptides_(.*)_allele_([0-9]+)_predictions.csv ]]
    uname="\${BASH_REMATCH[1]}_allele_\${BASH_REMATCH[2]}"
    echo \$uname

    concat_prediction.py \\
      -i $predictions \\
      -c ${chunk_size} \\
      -o predictions.buffer_\${uname}.tsv \\
      --pepmap "${PEPTIDE_MAP}" \\
      --allelemap "${ALLELE_MAP}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
