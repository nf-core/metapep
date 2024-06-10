process MERGE_PREDICTIONS_BUFFER {

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"


    input:
    path    predictions
    path    prediction_warnings

    output:
    path "predictions.buffer_*.tsv",            emit: ch_predictions_merged_buffer
    path "prediction_warnings.buffer_*.log",    emit: ch_prediction_warnings_merged_buffer
    path "versions.yml",                        emit: versions

    script:
    def chunk_size = params.prediction_chunk_size * params.pred_chunk_size_scaling
    """
    [[ ${predictions[0]} =~  peptides_(.*)_predictions.tsv ]];
    uname="\${BASH_REMATCH[1]}"
    echo \$uname

    concat_tsv.py -i $predictions -c $chunk_size -o predictions.buffer_\$uname.tsv
    sort -u $prediction_warnings > prediction_warnings.buffer_\$uname.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
