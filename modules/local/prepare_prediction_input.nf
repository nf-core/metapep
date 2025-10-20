process PREPARE_PREDICTION_INPUT {
  label 'process_single'
  tag "${meta.id}"   

  conda "bioconda::mhcgnomes=1.8.6 conda-forge::pandas=1.5.0"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.6--pyh7cba7a3_0' :
        'biocontainers/mhcgnomes:1.8.6--pyh7cba7a3_0' }"

  input:
  tuple val(meta), path(tsv)
  path  supported_alleles_json
  path  alleles_file

  output:
  tuple val(meta),
        path("${meta.id}_mhcflurry_input.csv"),
        path("${meta.id}_allele_supported.txt"),
        emit: flurry,  optional: true

  tuple val(meta),
        path("${meta.id}_mhcnuggets_input.tsv"),
        path("${meta.id}_allele_supported.txt"),
        emit: nuggets, optional: true

  //only for mhcnuggets-class-1 and mhcnuggets-class-2
  tuple val(meta),
        path("${meta.id}_idmap.csv"),
        emit: idmap,   optional: true

  path  "versions.yml", emit: versions

  script:
  """
  prepare_prediction_input.py \\
    --input ${tsv} \\
    --prefix ${meta.id} \\
    --pred_method ${params.pred_method} \\
    --allele_id "${meta.allele_id ?: ''}" \\
    --allele_name "${meta.allele_name ?: ''}" \\
    --mhc_class "${meta.mhc_class ?: params.mhc_class ?: ''}" \\
    --min_pep_len ${params.min_pep_len} \\
    --max_pep_len ${params.max_pep_len} \\
    --alleles ${alleles_file} \\
    --supported_alleles_json ${supported_alleles_json}
  """

  stub:
  def prefix = "${meta.id}"
  """
    touch ${prefix}_mhcflurry_input.csv
    touch ${prefix}_mhcnuggets_input.tsv
    touch ${prefix}_idmap.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        mhcgnomes: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcgnomes').version)")
    END_VERSIONS
  """
}