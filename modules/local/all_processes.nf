/*
 * Assign entity weights for input type 'assembly' and 'bins'
 */
process assign_nucl_entity_weights {
    publishDir "${params.outdir}/db_tables", mode: params.publish_dir_mode,
        saveAs: {filename -> "$filename" }

    input:
    val  microbiome_ids     from  ch_weights.microbiome_ids.collect().ifEmpty([])
    path weights_files      from  ch_weights.weights_paths.collect().ifEmpty([])

    output:
    path   "microbiomes_entities.nucl.tsv"    into   ch_nucl_microbiomes_entities  // entity_name, microbiome_id, entity_weight

    script:
    microbiome_ids = microbiome_ids.join(' ')
    """
    assign_entity_weights.py \
        --microbiome-ids $microbiome_ids \
        --weights-files $weights_files \
        --out microbiomes_entities.nucl.tsv
    """
}

/*
 * concat files and assign new, unique ids for all proteins (from different sources)
 */
process generate_protein_and_entity_ids {
    publishDir "${params.outdir}/db_tables", mode: params.publish_dir_mode,
        saveAs: {filename -> "$filename" }

    input:
    // Predicted Proteins
    path   predicted_proteins                  from       ch_pred_proteins.collect().ifEmpty([])
    val    predicted_proteins_microbiome_ids   from       ch_pred_proteins_microbiome_ids.collect().ifEmpty([])
    val    predicted_proteins_bin_basenames    from       ch_pred_proteins_bin_basename.collect().ifEmpty([])
    // Entrez Proteins
    path   entrez_proteins                     from       ch_entrez_proteins.ifEmpty([])
    path   entrez_entities_proteins            from       ch_entrez_entities_proteins.ifEmpty([])       //   protein_tmp_id (accessionVersion), entity_name (taxon_id)
    path   entrez_microbiomes_entities         from       ch_entrez_microbiomes_entities.ifEmpty([])    //   entity_name, microbiome_id, entity_weight
    // Bare Proteins
    path   bare_proteins                       from       ch_proteins_input.files.collect().ifEmpty([])
    path   bare_proteins_microbiome_ids        from       ch_proteins_input.ids.collect().ifEmpty([])

    output:
    path   "proteins.tsv.gz"                        into   ch_proteins
    path   "entities_proteins.tsv"                  into   ch_entities_proteins
    path   "entities.tsv"                           into   ch_entities
    path   "microbiomes_entities.no_weights.tsv"    into   ch_microbiomes_entities_noweights  // microbiome_id, entitiy_id  (no weights yet!)

    script:
    predicted_proteins_microbiome_ids = predicted_proteins_microbiome_ids.join(' ')
    predicted_proteins_bin_basenames  = predicted_proteins_bin_basenames.collect{ it ? it : "__ISASSEMBLY__" }.join(' ')
    """
    generate_protein_and_entity_ids.py \
        --predicted-proteins                  $predicted_proteins                  \
        --predicted-proteins-microbiome-ids   $predicted_proteins_microbiome_ids   \
        --predicted-proteins-bin-basenames    $predicted_proteins_bin_basenames    \
        --entrez-proteins                     "$entrez_proteins"                   \
        --entrez-entities-proteins            "$entrez_entities_proteins"          \
        --entrez-microbiomes-entities         "$entrez_microbiomes_entities"       \
        --bare-proteins                       $bare_proteins                       \
        --bare-proteins-microbiome-ids        $bare_proteins_microbiome_ids        \
        --out-proteins                        proteins.tsv.gz                      \
        --out-entities-proteins               entities_proteins.tsv                \
        --out-entities                        entities.tsv                         \
        --out-microbiomes-entities            microbiomes_entities.no_weights.tsv
    """
}

/*
 * Create microbiome_entities
 */
process finalize_microbiome_entities {
    publishDir "${params.outdir}/db_tables", mode: params.publish_dir_mode,
        saveAs: {filename -> "$filename" }

    input:
    path   entrez_microbiomes_entities        from       ch_entrez_microbiomes_entities.ifEmpty([])
    path   nucl_microbiomes_entities          from       ch_nucl_microbiomes_entities.ifEmpty([])
    path   microbiomes_entities_noweights     from       ch_microbiomes_entities_noweights
    path   entities                           from       ch_entities

    output:
    path   "microbiomes_entities.tsv"    into   ch_microbiomes_entities  // entity_id, microbiome_id, entity_weight

    script:

    """
    finalize_microbiome_entities.py \
        -eme $entrez_microbiomes_entities \
        -nme $nucl_microbiomes_entities \
        -menw $microbiomes_entities_noweights \
        -ent "$entities" \
        -o microbiomes_entities.tsv
    """
}

/*
 * Generate peptides
 */
process generate_peptides {
    publishDir "${params.outdir}", mode: params.publish_dir_mode,
        saveAs: {filename -> "db_tables/$filename" }

    input:
    file proteins from ch_proteins

    output:
    file "peptides.tsv.gz" into ch_peptides                // peptide_id, peptide_sequence
    file "proteins_peptides.tsv" into ch_proteins_peptides // protein_id, peptide_id, count
    //file "proteins_lengths.tsv"

    script:
    def min_pep_len = params.min_pep_len
    def max_pep_len = params.max_pep_len
    """
    generate_peptides.py -i $proteins \
                         -min $min_pep_len \
                         -max $max_pep_len \
                         -p "peptides.tsv.gz" \
                         -pp "proteins_peptides.tsv" \
                         -l "proteins_lengths.tsv"
    """
}

/*
 * Collect some numbers: proteins, peptides, unique peptides per conditon
 */
process collect_stats {
    publishDir "${params.outdir}", mode: params.publish_dir_mode,
        saveAs: {filename -> "db_tables/$filename" }

    input:
    path  peptides              from  ch_peptides
    path  proteins_peptides     from  ch_proteins_peptides
    path  entities_proteins     from  ch_entities_proteins
    path  microbiomes_entities  from  ch_microbiomes_entities
    path  conditions            from  ch_conditions

    output:
    file "stats.txt" into ch_stats

    script:
    """
    collect_stats.py --peptides "$peptides" \
                     --protein-peptide-occ "$proteins_peptides" \
                     --entities-proteins-occ "$entities_proteins" \
                     --microbiomes-entities-occ "$microbiomes_entities" \
                     --conditions "$conditions" \
                     --outfile stats.txt
    """
}

/*
 * Split prediction tasks (peptide, allele) into chunks of peptides that are to
 * be predicted against the same allele for parallel prediction
 */
process split_pred_tasks {
    input:
    path  peptides              from  ch_peptides
    path  proteins_peptides     from  ch_proteins_peptides
    path  entities_proteins     from  ch_entities_proteins
    path  microbiomes_entities  from  ch_microbiomes_entities
    path  conditions            from  ch_conditions
    path  conditions_alleles    from  ch_conditions_alleles
    path  alleles               from  ch_alleles
    // The tables are joined to map peptide -> protein -> microbiome -> condition -> allele
    // and thus to enumerate, which (peptide, allele) combinations have to be predicted.

    output:
    path "peptides_*.txt" into ch_epitope_prediction_chunks

    script:
    def pred_chunk_size       = params.pred_chunk_size
    def subsampling = params.sample_n ? "--sample_n ${params.sample_n}" : ""
    """
    gen_prediction_chunks.py --peptides "$peptides" \
                             --protein-peptide-occ "$proteins_peptides" \
                             --entities-proteins-occ "$entities_proteins" \
                             --microbiomes-entities-occ "$microbiomes_entities" \
                             --conditions "$conditions" \
                             --condition-allele-map "$conditions_alleles" \
                             --max-chunk-size $pred_chunk_size \
                             $subsampling \
                             --alleles "$alleles" \
                             --outdir .
    """
}

/*
 * Perform epitope prediction
 */
process predict_epitopes {
            // switch (params.pred_method) {
        // case "syfpeithi":
        //     pred_method_version = "1.0";
        //     break;
        // case "mhcflurry":
        //     pred_method_version = "1.4.3";
        //     break;
        // case "mhcnuggets-class-1":
        //     pred_method_version = "2.3.2";
        //     break;
        // case "mhcnuggets-class-2":
        //     pred_method_version = "2.3.2";
        //     break;
        // default:
        //     exit 1, "Epitope prediction method specified with --pred_method not recognized."
        // }
    input:
    path peptides from ch_epitope_prediction_chunks.flatten()

    output:
    path "*predictions.tsv" into ch_epitope_predictions
    path "*prediction_warnings.log" into ch_epitope_prediction_warnings

    script:
    def pred_method           = params.pred_method
    """

    # Extract allele name from file header
    allele_name="\$(head -n1 "$peptides" | fgrep '#' | cut -f2 -d'#')"
    allele_id="\$(head -n1 "$peptides" | fgrep '#' | cut -f3 -d'#')"

    out_basename="\$(basename "$peptides" .txt)"
    out_predictions="\$out_basename"_predictions.tsv
    out_warnings="\$out_basename"_prediction_warnings.log

    # Create output header
    echo "peptide_id	prediction_score	allele_id" >"\$out_predictions"

    # Process file
    # The --syfpeithi-norm flag enables score normalization when syfpeithi is
    # used and is ignored otherwise
    if ! epytope_predict.py --peptides "$peptides" \
                       --method "$pred_method" \
                       --method_version "$pred_method_version" \
		       --syfpeithi-norm \
                       "\$allele_name" \
                       2>stderr.log \
                       | tail -n +2 \
                       | cut -f 1,3 \
                       | sed -e "s/\$/	\$allele_id/" \
                       >>"\$out_basename"_predictions.tsv; then
        cat stderr.log >&2
        exit 1
    fi

    # Filter stderr for warnings and pass them on in the warnings channel
    fgrep WARNING stderr.log  | sort -u >"\$out_warnings" || :
    """
}

/*
 * Merge prediction results from peptide chunks into one prediction result
 */
 // gather chunks of predictions and merge them already to avoid too many input files for `merge_predictions` process
 // (causing "sbatch: error: Batch job submission failed: Pathname of a file, directory or other parameter too long")
 // sort and buffer to ensure resume will work (inefficient, since this causes waiting for all predictions)
ch_epitope_predictions_buffered = ch_epitope_predictions.toSortedList().flatten().buffer(size: 1000, remainder: true)
ch_epitope_prediction_warnings_buffered = ch_epitope_prediction_warnings.toSortedList().flatten().buffer(size: 1000, remainder: true)

process merge_predictions_buffer {

    input:
    path predictions from ch_epitope_predictions_buffered
    path prediction_warnings from ch_epitope_prediction_warnings_buffered

    output:
    path "predictions.buffer_*.tsv" into ch_predictions_merged_buffer
    path "prediction_warnings.buffer_*.log" into ch_prediction_warnings_merged_buffer

    script:
    def single = predictions instanceof Path ? 1 : predictions.size()
    def merge = (single == 1) ? 'cat' : 'csvtk concat -t'
    """
    [[ ${predictions[0]} =~  peptides_(.*)_predictions.tsv ]];
    uname="\${BASH_REMATCH[1]}"
    echo \$uname

    $merge $predictions > predictions.buffer_\$uname.tsv
    sort -u $prediction_warnings > prediction_warnings.buffer_\$uname.log
    """
}

process merge_predictions {
    publishDir "${params.outdir}", mode: params.publish_dir_mode,
        saveAs: {filename -> filename.endsWith(".log") ? "logs/$filename" : "db_tables/$filename"}

    input:
    path predictions from ch_predictions_merged_buffer.collect()
    path prediction_warnings from ch_prediction_warnings_merged_buffer.collect()

    output:
    path "predictions.tsv.gz" into ch_predictions
    path "prediction_warnings.log"

    script:
    def single = predictions instanceof Path ? 1 : predictions.size()
    def merge = (single == 1) ? 'cat' : 'csvtk concat -t'
    """
    $merge $predictions | gzip > predictions.tsv.gz
    sort -u $prediction_warnings > prediction_warnings.log
    """
}

/*
 * Generate figures
 */
process prepare_score_distribution {
    publishDir "${params.outdir}/figures/prediction_scores", mode: params.publish_dir_mode

    input:
    file predictions from ch_predictions
    file proteins_peptides from ch_proteins_peptides
    file entities_proteins from ch_entities_proteins
    file microbiomes_entities from ch_microbiomes_entities
    file conditions from  ch_conditions
    file conditions_alleles from  ch_conditions_alleles
    file alleles from ch_alleles

    output:
    file "prediction_scores.allele_*.tsv" into ch_prep_prediction_scores

    script:
    """
    prepare_score_distribution.py --predictions "$predictions" \
                            --protein-peptide-occ "$proteins_peptides" \
                            --entities-proteins-occ "$entities_proteins" \
                            --microbiomes-entities-occ "$microbiomes_entities" \
                            --conditions "$conditions" \
                            --condition-allele-map "$conditions_alleles" \
                            --alleles "$alleles" \
                            --outdir .
    """
}

process plot_score_distribution {
    publishDir "${params.outdir}/figures", mode: params.publish_dir_mode

    input:
    file prep_scores from ch_prep_prediction_scores.flatten()
    file alleles from ch_alleles
    file conditions from ch_conditions

    output:
    file "prediction_score_distribution.*.pdf"

    script:
    """
    [[ ${prep_scores} =~ prediction_scores.allele_(.*).tsv ]];
    allele_id="\${BASH_REMATCH[1]}"
    echo \$allele_id

    plot_score_distribution.R --scores $prep_scores \
                                   --alleles $alleles \
                                   --conditions $conditions \
                                   --allele_id \$allele_id \
                                   --method ${params.pred_method}
    """

}

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