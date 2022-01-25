process MERGE_PREDICTIONS_BUFFER {
    label 'cache_lenient'
    label 'process_medium_memory'

    conda (params.enable_conda ? "bioconda::csvtk=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0' :
        'quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0' }"

    input:
    path    predictions
    path    prediction_warnings

    output:
    path "predictions.buffer_*.tsv",            emit: ch_predictions_merged_buffer
    path "prediction_warnings.buffer_*.log",    emit: ch_prediction_warnings_merged_buffer
    path "versions.yml",                        emit: versions

    script:
    def single = predictions instanceof Path ? 1 : predictions.size()
    def merge = (single == 1) ? 'cat' : 'csvtk concat -t'
    """
    [[ ${predictions[0]} =~  peptides_(.*)_predictions.tsv ]];
    uname="\${BASH_REMATCH[1]}"
    echo \$uname

    $merge $predictions > predictions.buffer_\$uname.tsv
    sort -u $prediction_warnings > prediction_warnings.buffer_\$uname.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
