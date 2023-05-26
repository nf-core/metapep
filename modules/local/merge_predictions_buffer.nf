process MERGE_PREDICTIONS_BUFFER {
    label 'cache_lenient'
    label 'process_medium_memory'

    // TODO generate extra biocontainer with only specifying pandas version (currently mulled container taken from "bedtools=2.23.0,pandas=1.5.2")
    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d19e2715c83e4582e3f1fb0a2e473abde8ca636e:fc171b36fc2e2a38a259a1c82a139b59d94c968b-0' :
        'biocontainers/mulled-v2-d19e2715c83e4582e3f1fb0a2e473abde8ca636e:fc171b36fc2e2a38a259a1c82a139b59d94c968b-0' }"

    input:
    path    predictions
    path    prediction_warnings

    output:
    path "predictions.buffer_*.tsv",            emit: ch_predictions_merged_buffer
    path "prediction_warnings.buffer_*.log",    emit: ch_prediction_warnings_merged_buffer
    path "versions.yml",                        emit: versions

    script:
    def chunk_size = params.ds_prep_chunk_size
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
