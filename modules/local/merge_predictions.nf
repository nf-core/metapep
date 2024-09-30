process MERGE_PREDICTIONS {
    label "process_long"

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"


    input:
    path predictions
    path prediction_warnings

    output:
    path "predictions.tsv.gz",          emit: ch_predictions
    path "prediction_warnings.log",     emit: ch_prediction_warnings
    path "versions.yml",                emit: versions

    script:
    def chunk_size = params.prediction_chunk_size * params.pred_chunk_size_scaling
    """
    concat_tsv.py -i $predictions -c $chunk_size -o predictions.tsv.gz
    sort -u $prediction_warnings > prediction_warnings.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
