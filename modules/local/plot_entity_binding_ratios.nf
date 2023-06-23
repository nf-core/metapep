process PLOT_ENTITY_BINDING_RATIOS {
    label 'cache_lenient'
    label 'process_medium_memory'

    conda "conda-forge::r-ggplot2=3.4.2 conda-forge::r-data.table=1.14.8 conda-forge::r-dplyr=1.1.2 conda-forge::r-stringr=1.5.0 conda-forge::r-ggpubr=0.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ea91ce9d3f5052e6d1ab53975ea9774ed365996c:5a35843d615168187ed1757a15ea488ee0a0d634-0' :
        'biocontainers/mulled-v2-ea91ce9d3f5052e6d1ab53975ea9774ed365996c:5a35843d615168187ed1757a15ea488ee0a0d634-0' }"

    publishDir "${params.outdir}/figures", mode: params.publish_dir_mode

    input:
    path prep_entity_binding_ratios
    path alleles

    output:
    path "entity_binding_ratios.*.pdf",     emit:   ch_plot_entity_binding_ratios
    path "versions.yml",                    emit:   versions

    script:
    def handle_pvalue = params.show_pvalue ? "-t TRUE" : "-t FALSE"
    """
    [[ ${prep_entity_binding_ratios} =~ entity_binding_ratios.allele_(.*).tsv ]];
    allele_id="\${BASH_REMATCH[1]}"
    echo \$allele_id

    plot_entity_binding_ratios.R \\
        -r $prep_entity_binding_ratios \\
        -a $alleles \\
        $handle_pvalue \\
        \$allele_id

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
