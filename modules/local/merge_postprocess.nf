process MERGE_POSTPROCESS {
  label 'process_short'
  conda "conda-forge::pandas=1.5.2"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
      'biocontainers/pandas:1.5.2' }"

  input:
  path tsvs

  output:
  path "predictions.tsv.gz", emit: predictions
  path "versions.yml",      emit: versions

  script:
  """
  merge_postprocess.py \\
      --inputs ${tsvs.join(' ')} \\
      --out predictions.tsv.gz

  # versions.yml erstellen
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    python: \$(python --version | sed 's/Python //g')
    pandas: \$(python -c "import pandas as pd; print(pd.__version__)")
  END_VERSIONS  
  """
}
