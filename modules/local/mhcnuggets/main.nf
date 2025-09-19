process MHCNUGGETS {
  label 'process_single'
  tag "${meta.id}"

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/mhcnuggets:2.4.0--pyh7cba7a3_0' :
      'quay.io/biocontainers/mhcnuggets:2.4.0--pyh7cba7a3_0' }"

  input:
  tuple val(meta), path(tsv)

  output:
  tuple val(meta), path("*_predicted_mhcnuggets*.csv"), emit: predicted
  path  "versions.yml",                                   emit: versions

  // WICHTIG: nur Template
  script:
  template "mhcnuggets.py"

  stub:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo "peptide,allele,ic50,rank" > ${prefix}_predicted_mhcnuggets.csv
  echo "AAAA,H-2-Db,50,0.1"       >> ${prefix}_predicted_mhcnuggets.csv
  echo -e "MHCNUGGETS:\\n  pandas: 1.5.2" > versions.yml
  """
}
