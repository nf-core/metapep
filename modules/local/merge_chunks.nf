process MERGE_CHUNKS {
    label "process_long"

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path predictions
    path peptide_map
    path allele_map

    output:
    path "predictions.tsv.gz"     , emit: ch_predictions
    path "versions.yml"           , emit: versions

    script:
    def chunk_size = params.prediction_chunk_size * params.pred_chunk_size_scaling
    """

    concat_prediction.py \\
        -i $predictions \\
        -o predictions.tsv.gz \\
        -c $chunk_size \\
        --peptides "${peptide_map}" \\
        --alleles "${allele_map}"

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """
}
