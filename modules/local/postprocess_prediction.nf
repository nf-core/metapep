process POSTPROCESS_PREDICTIONS {
  label 'process_short'
  tag   "${meta.id}"

  conda "conda-forge::pandas=1.5.2"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
      'biocontainers/pandas:1.5.2' }"

  input:
  tuple val(meta), val(method), path(input_csv), path(predicted_csv)

  output:
  path "${meta.id}_${method}_reduced.tsv", emit: reduced
  path "versions.yml",                     emit: versions

  script:
  def outFile = "${meta.id}_${method}_reduced.tsv"
  """
  postprocess_prediction.py \\
      --method ${method} \\
      --input  "${input_csv}" \\
      --pred   "${predicted_csv}" \\
      --out    "${outFile}" \\
      
  
  # versions.yml wie in den alten Modulen schreiben
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      python: \$(python --version | sed 's/Python //g')
      pandas: \$(python -c "import pandas as pd; print(pd.__version__)")
      tool: ${method}
  END_VERSIONS
  """

}
