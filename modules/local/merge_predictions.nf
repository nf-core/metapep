process MERGE_PREDICTIONS {
    label "process_high_memory"
    label 'cache_lenient'

    conda "bioconda::csvtk=0.23.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0' :
        'quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0' }"

    input:
    path predictions
    path prediction_warnings

    output:
    path "predictions.tsv.gz",          emit: ch_predictions
    path "prediction_warnings.log",     emit: ch_prediction_warnings
    path "versions.yml",                emit: versions

    script:
    def single = predictions instanceof Path ? 1 : predictions.size()
    def merge = (single == 1) ? 'cat' : 'csvtk concat -t'
    """
    $merge $predictions | gzip > predictions.tsv.gz
    sort -u $prediction_warnings > prediction_warnings.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
