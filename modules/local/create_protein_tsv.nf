process CREATE_PROTEIN_TSV {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::biopython=1.78"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(protein_fasta)

    output:
    tuple val(meta), path("proteins.pred_${meta.id}*.tsv.gz"), emit: ch_pred_proteins     // Emit protein tsv
    path    "versions.yml"                                   , emit: versions

    script:
    name = meta.bin_basename ? "${meta.id}.${meta.bin_basename}" : "${meta.id}"
    """
    echo -e "protein_tmp_id\tprotein_sequence" > proteins.pred_${name}.tsv
    fasta_to_tsv.py --remove-asterisk --input ${protein_fasta} >> proteins.pred_${name}.tsv
    gzip proteins.pred_${name}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
