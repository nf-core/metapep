process PLOT_SCORE_DISTRIBUTION {
    label 'cache_lenient'
    label 'process_medium_memory'

    conda "bioconda::bioconductor-alphabeta=1.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-alphabeta:1.8.0--r41hdfd78af_0' :
        'biocontainers/bioconductor-alphabeta:1.8.0--r41hdfd78af_0' }"

    input:
    path prep_scores
    path alleles
    path conditions

    output:
    path "prediction_score_distribution.*.pdf",     emit:   ch_plot_score_distribution
    path "versions.yml",                            emit:   versions

    script:
    def syfpeithi_threshold = params.syfpeithi_score_threshold
    def mhcfn_threshold = params.mhcflurry_mhcnuggets_score_threshold
    """
    [[ ${prep_scores} =~ prediction_scores.allele_(.*).tsv ]];
    allele_id="\${BASH_REMATCH[1]}"
    echo \$allele_id

    plot_score_distribution.R \\
        $prep_scores \\
        $alleles \\
        $conditions \\
        \$allele_id \\
        ${params.pred_method} \\
        $syfpeithi_threshold \\
        $mhcfn_threshold


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        stringr: \$(Rscript -e "library(stringr); cat(as.character(packageVersion('stringr')))")
    END_VERSIONS
    """

}
