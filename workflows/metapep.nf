/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PRODIGAL               } from '../modules/nf-core/prodigal/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metapep_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SHOW_SUPPORTED_MODELS             } from '../modules/local/show_supported_models'
include { DOWNLOAD_PROTEINS                 } from '../modules/local/download_proteins'
include { CREATE_PROTEIN_TSV                } from '../modules/local/create_protein_tsv'
include { ASSIGN_NUCL_ENTITY_WEIGHTS        } from '../modules/local/assign_nucl_entity_weights'
include { GENERATE_PROTEIN_AND_ENTITY_IDS   } from '../modules/local/generate_protein_and_entity_ids'
include { FINALIZE_MICROBIOME_ENTITIES      } from '../modules/local/finalize_microbiome_entities'
include { GENERATE_PEPTIDES                 } from '../modules/local/generate_peptides'
include { COLLECT_STATS                     } from '../modules/local/collect_stats'
include { SPLIT_PRED_TASKS                  } from '../modules/local/split_pred_tasks'
include { MERGE_CHUNKS_BUFFER               } from '../modules/local/merge_chunks_buffer'
include { MERGE_CHUNKS                      } from '../modules/local/merge_chunks'
include { PREPARE_SCORE_DISTRIBUTION        } from '../modules/local/prepare_score_distribution'
include { PLOT_SCORE_DISTRIBUTION           } from '../modules/local/plot_score_distribution'
include { PREPARE_ENTITY_BINDING_RATIOS     } from '../modules/local/prepare_entity_binding_ratios'
include { PLOT_ENTITY_BINDING_RATIOS        } from '../modules/local/plot_entity_binding_ratios'

include { PROCESS_INPUT                     } from '../subworkflows/local/process_input'

// Subworkflow (implemented from nf-core/epitopeprediction)
include { MHC_BINDING_PREDICTION            } from '../subworkflows/local/mhc_binding_prediction'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAPEP {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Load supported alleles file
    supported_alleles_json = file("$projectDir/assets/supported_alleles.json", checkIfExists: true)
    netmhc_software_meta   = file("$projectDir/assets/netmhc_software_meta.json", checkIfExists: true)

    //
    // SUBBRANCH: Show supported alleles for all prediction methods
    //
    if (params.show_supported_models) {
        ch_supported_alleles_and_peptide_lengths = Channel.fromPath(
            supported_alleles_json
        )
        SHOW_SUPPORTED_MODELS(ch_supported_alleles_and_peptide_lengths)
        ch_versions = ch_versions.mix(SHOW_SUPPORTED_MODELS.out.versions)
    } else {

        //
        // SUBWORKFLOW: Read in samplesheet, validate and create/populate datamodel based on input
        //
        PROCESS_INPUT (
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(PROCESS_INPUT.out.versions)

        //
        // MODULE: Download proteins from entrez
        //
        DOWNLOAD_PROTEINS (
            PROCESS_INPUT.out.ch_taxa_input.map { meta, file -> meta.id }.collect(),
            PROCESS_INPUT.out.ch_taxa_input.map { meta, file -> file }.collect()
        )
        ch_versions = ch_versions.mix(DOWNLOAD_PROTEINS.out.versions)

        //
        // MODULE: Predict proteins from nucleotides
        //
        PRODIGAL(
            PROCESS_INPUT.out.ch_nucl_input,
            "gff"
        )
        ch_versions = ch_versions.mix(PRODIGAL.out.versions)

        CREATE_PROTEIN_TSV (
            PRODIGAL.out.amino_acid_fasta
        )
        ch_versions = ch_versions.mix(CREATE_PROTEIN_TSV.out.versions)

        //
        // MODULE: Assign entity weights (Nucleotide Input)
        //
        ASSIGN_NUCL_ENTITY_WEIGHTS (
            PROCESS_INPUT.out.ch_weights.map { meta, file -> meta.id }.collect().ifEmpty([]),
            PROCESS_INPUT.out.ch_weights.map { meta, file -> file }.collect().ifEmpty([])
        )
        ch_versions = ch_versions.mix(ASSIGN_NUCL_ENTITY_WEIGHTS.out.versions)

        //
        // MODULE: Generate protein and entity ids
        //
        // concat files and assign new, unique ids for all proteins (from different sources)
        // Sort predicted protein input for GENERATE_PROTEIN_AND_ENTITY_IDS to ensure deterministic id assignments
        ch_pred_proteins_sorted = CREATE_PROTEIN_TSV.out.ch_pred_proteins.toSortedList( { a, b -> a[0].id <=> b[0].id } ).flatMap()

        GENERATE_PROTEIN_AND_ENTITY_IDS (
            PROCESS_INPUT.out.ch_microbiomes,
            ch_pred_proteins_sorted.collect { meta, file -> file }.ifEmpty([]),
            ch_pred_proteins_sorted.collect { meta, file -> meta }.ifEmpty([]),
            DOWNLOAD_PROTEINS.out.ch_entrez_proteins.ifEmpty([]),
            DOWNLOAD_PROTEINS.out.ch_entrez_entities_proteins.ifEmpty([]),
            DOWNLOAD_PROTEINS.out.ch_entrez_microbiomes_entities.ifEmpty([])
        )
        ch_versions = ch_versions.mix(GENERATE_PROTEIN_AND_ENTITY_IDS.out.versions)

        //
        // MODULE: Create microbiome_entities
        //
        FINALIZE_MICROBIOME_ENTITIES (
            DOWNLOAD_PROTEINS.out.ch_entrez_microbiomes_entities.ifEmpty([]),
            ASSIGN_NUCL_ENTITY_WEIGHTS.out.ch_nucl_microbiomes_entities.ifEmpty([]),
            GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_microbiomes_entities_noweights,
            GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities
        )
        ch_versions = ch_versions.mix(FINALIZE_MICROBIOME_ENTITIES.out.versions)

        //
        // MODULE: Generate peptides
        //
        GENERATE_PEPTIDES (
            GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_proteins,
            PROCESS_INPUT.out.peptide_lengths
        )
        ch_versions = ch_versions.mix(GENERATE_PEPTIDES.out.versions)

        //
        // MODULE: Split prediction tasks into chunks
        //

        // Split prediction tasks (peptide, allele) into chunks of peptides that are to
        // be predicted against the same allele for parallel prediction
        SPLIT_PRED_TASKS (
            GENERATE_PEPTIDES.out.ch_peptides,
            GENERATE_PEPTIDES.out.ch_proteins_peptides,
            GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities_proteins,
            FINALIZE_MICROBIOME_ENTITIES.out.ch_microbiomes_entities,
            PROCESS_INPUT.out.ch_conditions,
            PROCESS_INPUT.out.ch_conditions_alleles,
            PROCESS_INPUT.out.ch_alleles
        )
        ch_versions = ch_versions.mix(SPLIT_PRED_TASKS.out.versions)

        //
        // SUBWORKFLOW: MHC binding prediction
        //
        // Build tuples [meta, TSV] by parsing the header of the TXT, but use the .tsv as file
        SPLIT_PRED_TASKS.out.ch_epitope_prediction_chunks
            .flatten()
            .map { txt -> tuple(txt.baseName, txt) }
            .join(
                SPLIT_PRED_TASKS.out.ch_epitope_prediction_chunks_tsv
                    .flatten()
                    .map { tsv -> tuple(tsv.baseName, tsv) },
                by: 0
            )
            .map { base, txt, tsv ->
                // Parse header inline
                def lines = txt.text.readLines()
                def header = lines.isEmpty() ? '' : lines[0]

                if (!header.startsWith('#')) {
                    return tuple(base, txt, tsv, null)
                }

                def parts = header.substring(1).split('#', 3)
                def allele_name = parts[0]?.trim()
                def allele_id = parts.size() > 1 ? parts[1]?.trim() : ''

                if (!allele_name || allele_name.equalsIgnoreCase('NA') || allele_name.equalsIgnoreCase('null')) {
                    return tuple(base, txt, tsv, null)
                }

                def header_info = [allele_name: allele_name, allele_id: allele_id]
                tuple(base, txt, tsv, header_info)
            }
            .filter { base, txt, tsv, header_info -> header_info != null }
            .map { base, txt, tsv, header_info ->
                def mhc_class = params.pred_method in ['mhcnuggets', 'mhcflurry', 'netmhcpan'] ? 'I' : 'II'

                tuple(
                    [
                        id          : "${base}_allele_${header_info.allele_id}",
                        file_id     : base,
                        allele_id   : header_info.allele_id,
                        allele_name : header_info.allele_name,
                        alleles     : header_info.allele_name,
                        mhc_class   : mhc_class
                    ],
                    tsv
                )
            }
            .set { ch_chunks_with_meta }

        MHC_BINDING_PREDICTION(
            ch_chunks_with_meta,
            params.pred_method,
            supported_alleles_json,
            netmhc_software_meta
        )
        ch_versions = ch_versions.mix(MHC_BINDING_PREDICTION.out.versions)

        //
        // MODULE: Merge prediction results
        //

        // Count all generated files for prediction for buffering decision
        // MHC_BINDING_PREDICTION.out.predicted is tuple val(meta), path(csv)
        // Extract just the files for counting and merging
        // Generates a channel with tuples: [count, prediction_file]
        // the channel is branched into buffer or unbuffer depending on the file count
        // therefore one of both channels will be empty
        SPLIT_PRED_TASKS.out.ch_epitope_prediction_chunks.flatten().count()
            .combine(MHC_BINDING_PREDICTION.out.predicted.map { meta, file -> file }.toSortedList().flatten())
            .branch { count, predictions_file ->
                buffer: count > params.pred_buffer_files
                    return predictions_file
                unbuffered: count <= params.pred_buffer_files
                    return predictions_file
            }.set { ch_pred_merge_input }

        // Process is only used when files exceed the buffer files parameter (default 1000) -> May generate issues for slurm if larger
        MERGE_CHUNKS_BUFFER (
            ch_pred_merge_input.buffer.buffer(size: params.pred_buffer_files, remainder: true),
            GENERATE_PEPTIDES.out.ch_peptides.first(),
            PROCESS_INPUT.out.ch_alleles.first()
        )
        ch_versions = ch_versions.mix(MERGE_CHUNKS_BUFFER.out.versions)

        // Mix the output of the merge predictions buffer channel and merge predictions channel (one of them will be empty)
        ch_merge_predictions_input = MERGE_CHUNKS_BUFFER.out.ch_predictions_merged_buffer.mix(ch_pred_merge_input.unbuffered)

        MERGE_CHUNKS (
            ch_merge_predictions_input.collect(sort: { it.baseName }),
            GENERATE_PEPTIDES.out.ch_peptides.first(),
            PROCESS_INPUT.out.ch_alleles.first()
        )
        ch_versions = ch_versions.mix(MERGE_CHUNKS.out.versions)

        //
        // MODULE: Collect stats
        //
        // Collects proteins, peptides, unique peptides per condition
        COLLECT_STATS (
            GENERATE_PEPTIDES.out.ch_proteins_peptides,
            GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities_proteins,
            FINALIZE_MICROBIOME_ENTITIES.out.ch_microbiomes_entities,
            PROCESS_INPUT.out.ch_conditions,
            PROCESS_INPUT.out.ch_alleles,
            PROCESS_INPUT.out.ch_conditions_alleles,
            MERGE_CHUNKS.out.ch_predictions
        )
        ch_versions = ch_versions.mix(COLLECT_STATS.out.versions)

        //
        // MODULE: Plot score distributions
        //
        PREPARE_SCORE_DISTRIBUTION (
            MERGE_CHUNKS.out.ch_predictions,
            GENERATE_PEPTIDES.out.ch_proteins_peptides,
            GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities_proteins,
            FINALIZE_MICROBIOME_ENTITIES.out.ch_microbiomes_entities,
            PROCESS_INPUT.out.ch_conditions,
            PROCESS_INPUT.out.ch_conditions_alleles,
            PROCESS_INPUT.out.ch_alleles
        )
        ch_versions = ch_versions.mix(PREPARE_SCORE_DISTRIBUTION.out.versions)

        PLOT_SCORE_DISTRIBUTION (
            PREPARE_SCORE_DISTRIBUTION.out.ch_prep_prediction_scores.flatten(),
            PROCESS_INPUT.out.ch_alleles,
            PROCESS_INPUT.out.ch_conditions
        )
        ch_versions = ch_versions.mix(PLOT_SCORE_DISTRIBUTION.out.versions)

        //
        // MODULE: Plot entity binding ratios
        //
        PREPARE_ENTITY_BINDING_RATIOS (
            MERGE_CHUNKS.out.ch_predictions,
            GENERATE_PEPTIDES.out.ch_proteins_peptides,
            GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities_proteins,
            FINALIZE_MICROBIOME_ENTITIES.out.ch_microbiomes_entities,
            PROCESS_INPUT.out.ch_conditions,
            PROCESS_INPUT.out.ch_conditions_alleles,
            PROCESS_INPUT.out.ch_alleles
        )
        ch_versions = ch_versions.mix(PREPARE_ENTITY_BINDING_RATIOS.out.versions)

        PLOT_ENTITY_BINDING_RATIOS (
            PREPARE_ENTITY_BINDING_RATIOS.out.ch_prep_entity_binding_ratios.flatten(),
            PROCESS_INPUT.out.ch_alleles
        )
        ch_versions = ch_versions.mix(PLOT_ENTITY_BINDING_RATIOS.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'metapep_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.fromPath("$projectDir/assets/nf-core-metapep_logo_light.png")

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    if (!params.show_supported_models) {
    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_STATS.out.ch_stats)
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    multiqc_report = MULTIQC.out.report.toList()

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
