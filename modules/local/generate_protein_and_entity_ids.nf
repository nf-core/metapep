process GENERATE_PROTEIN_AND_ENTITY_IDS {
    label 'process_low'

    conda "conda-forge::pandas=1.5.2 conda-forge::biopython=1.79 conda-forge::numpy=1.23.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0' :
        'biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0' }"

    input:
    path(microbiomes)
    path(predicted_proteins)
    val(predicted_proteins_meta)
    path(entrez_proteins)
    path(entrez_entities_proteins)
    path(entrez_microbiomes_entities)

    output:
    path   "proteins.tsv.gz"                        , emit:   ch_proteins
    path   "entities_proteins.tsv"                  , emit:   ch_entities_proteins
    path   "entities.tsv"                           , emit:   ch_entities
    path   "microbiomes_entities.no_weights.tsv"    , emit:   ch_microbiomes_entities_noweights  // microbiome_id, entitiy_id  (no weights yet!)
    path   "versions.yml"                           , emit:   versions

    script:
        predicted_proteins_microbiome_ids   = predicted_proteins_meta.collect { meta -> meta.id }.join(' ')
        predicted_proteins_bin_basenames    = predicted_proteins_meta.collect { meta -> meta.bin_basename ?: "__ISASSEMBLY__" }.join(' ')

    """
    generate_protein_and_entity_ids.py \
        --microbiomes                         $microbiomes                         \\
        --predicted-proteins                  $predicted_proteins                  \\
        --predicted-proteins-microbiome-ids   $predicted_proteins_microbiome_ids   \\
        --predicted-proteins-bin-basenames    $predicted_proteins_bin_basenames    \\
        --entrez-proteins                     "$entrez_proteins"                   \\
        --entrez-entities-proteins            "$entrez_entities_proteins"          \\
        --entrez-microbiomes-entities         "$entrez_microbiomes_entities"       \\
        --out-proteins                        proteins.tsv.gz                      \\
        --out-entities-proteins               entities_proteins.tsv                \\
        --out-entities                        entities.tsv                         \\
        --out-microbiomes-entities            microbiomes_entities.no_weights.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
