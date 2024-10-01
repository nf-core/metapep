process GENERATE_PEPTIDES {
    label 'process_long'
    label 'process_medium_memory'
    label 'cache_lenient'

    conda "conda-forge::pandas=1.5.2 conda-forge::biopython=1.79 conda-forge::numpy=1.23.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0' :
        'biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0' }"

    input:
    path(proteins)
    val(peptide_lengths)

    output:
    path "peptides.tsv.gz"      , emit: ch_peptides               // peptide_id, peptide_sequence
    path "proteins_peptides.tsv", emit: ch_proteins_peptides      // protein_id, peptide_id, count
    path "versions.yml"         , emit: versions
    //file "proteins_lengths.tsv"

    script:
    def mem_log_level = params.memory_usage_log_deep ? "--mem_log_level_deep" : ""
    """
    generate_peptides.py -i $proteins \\
                        -p "peptides.tsv.gz" \\
                        -pp "proteins_peptides.tsv" \\
                        -l "proteins_lengths.tsv" \\
                        $mem_log_level \\
                        -pll ${peptide_lengths.join(" ")}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
