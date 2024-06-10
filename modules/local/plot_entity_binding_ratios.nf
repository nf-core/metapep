process PLOT_ENTITY_BINDING_RATIOS {
    label 'process_medium_memory'

    conda "conda-forge::r-ggplot2=3.4.2 conda-forge::r-data.table=1.14.8 conda-forge::r-dplyr=1.1.2 conda-forge::r-stringr=1.5.0 conda-forge::r-ggpubr=0.6.0 conda-forge::r-optparse=1.7.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0be74e7b0c2e289bc8098b1491baf4f181012b1c:a1635746bc2c13635cbea8c29bd5a2837bdd7cd5-0' :
        'biocontainers/mulled-v2-0be74e7b0c2e289bc8098b1491baf4f181012b1c:a1635746bc2c13635cbea8c29bd5a2837bdd7cd5-0' }"

    publishDir "${params.outdir}/figures", mode: params.publish_dir_mode

    input:
    path prep_entity_binding_ratios
    path alleles

    output:
    path "entity_binding_ratios.*.pdf",     emit:   ch_plot_entity_binding_ratios
    path "versions.yml",                    emit:   versions

    script:
    def hide_pvalue = params.hide_pvalue ? "TRUE" : "FALSE"
    """
    [[ ${prep_entity_binding_ratios} =~ entity_binding_ratios.allele_(.*).tsv ]];
    allele_id="\${BASH_REMATCH[1]}"
    echo \$allele_id

    plot_entity_binding_ratios.R \\
        -r $prep_entity_binding_ratios \\
        -a $alleles \\
        -d $hide_pvalue \\
        -i \$allele_id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        stringr: \$(Rscript -e "library(stringr); cat(as.character(packageVersion('stringr')))")
        ggpubr: \$(Rscript -e "library(ggpubr); cat(as.character(packageVersion('ggpubr')))")
        optparse: \$(Rscript -e "library(optparse); cat(as.character(packageVersion('optparse')))")
    END_VERSIONS
    """
}
