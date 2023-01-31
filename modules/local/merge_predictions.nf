process MERGE_PREDICTIONS {
    label "process_high_memory"
    label 'cache_lenient'

    // TODO generate extra biocontainer with only specifying pandas version (currently mulled container taken from "bedtools=2.23.0,pandas=1.5.2")
    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d19e2715c83e4582e3f1fb0a2e473abde8ca636e:fc171b36fc2e2a38a259a1c82a139b59d94c968b-0' :
        'quay.io/biocontainers/mulled-v2-d19e2715c83e4582e3f1fb0a2e473abde8ca636e:fc171b36fc2e2a38a259a1c82a139b59d94c968b-0' }"

    input:
    path predictions
    path prediction_warnings

    output:
    path "predictions.tsv.gz",          emit: ch_predictions
    path "prediction_warnings.log",     emit: ch_prediction_warnings
    path "versions.yml",                emit: versions

    script:
    """
    concat_tsv.py -i $predictions -c 100000 -o predictions.tsv.gz
    sort -u $prediction_warnings > prediction_warnings.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
