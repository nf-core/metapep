process prepare_entity_binding_ratios {
    publishDir "${params.outdir}/figures/entity_binding_ratios", mode: params.publish_dir_mode

    input:
    file predictions from ch_predictions
    file proteins_peptides from ch_proteins_peptides
    file entities_proteins from ch_entities_proteins
    file microbiomes_entities from ch_microbiomes_entities
    file conditions from  ch_conditions
    file conditions_alleles from  ch_conditions_alleles
    file alleles from ch_alleles

    output:
    file "entity_binding_ratios.allele_*.tsv" into ch_prep_entity_binding_ratios

    script:
    """
    prepare_entity_binding_ratios.py --predictions "$predictions" \
                            --protein-peptide-occ "$proteins_peptides" \
                            --entities-proteins-occ "$entities_proteins" \
                            --microbiomes-entities-occ "$microbiomes_entities" \
                            --conditions "$conditions" \
                            --condition-allele-map "$conditions_alleles" \
                            --alleles "$alleles" \
                            --method ${params.pred_method} \
                            --outdir .
    """
}

process plot_entity_binding_ratios {
    publishDir "${params.outdir}/figures", mode: params.publish_dir_mode

    input:
    file prep_entity_binding_ratios from ch_prep_entity_binding_ratios.flatten()
    file alleles from ch_alleles

    output:
    file "entity_binding_ratios.*.pdf"

    script:
    """
    [[ ${prep_entity_binding_ratios} =~ entity_binding_ratios.allele_(.*).tsv ]];
    allele_id="\${BASH_REMATCH[1]}"
    echo \$allele_id

    plot_entity_binding_ratios.R --binding-rates $prep_entity_binding_ratios \
                                   --alleles $alleles \
                                   --allele_id \$allele_id
    """
}