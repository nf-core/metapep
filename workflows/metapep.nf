/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetapep.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_PROTEINS                 } from '../modules/local/download_proteins'
include { CREATE_PROTEIN_TSV                } from '../modules/local/create_protein_tsv'
include { ASSIGN_NUCL_ENTITY_WEIGHTS        } from '../modules/local/assign_nucl_entity_weights'
include { GENERATE_PROTEIN_AND_ENTITY_IDS   } from '../modules/local/generate_protein_and_entity_ids'
include { FINALIZE_CONDITION_ENTITIES       } from '../modules/local/finalize_condition_entities'
include { FRED2_GENERATEPEPTIDES            } from '../modules/local/fred2_generatepeptides.nf'
include { SPLIT_PEPTIDES                    } from '../modules/local/split_peptides.nf'
include { GENERATE_PEPTIDES                 } from '../modules/local/generate_peptides'
include { COLLECT_STATS                     } from '../modules/local/collect_stats'
include { SPLIT_PRED_TASKS                  } from '../modules/local/split_pred_tasks'
include { PREDICT_EPITOPES                  } from '../modules/local/predict_epitopes'
include { MERGE_PREDICTIONS_BUFFER          } from '../modules/local/merge_predictions_buffer'
include { MERGE_PREDICTIONS                 } from '../modules/local/merge_predictions'
include { PREPARE_SCORE_DISTRIBUTION        } from '../modules/local/prepare_score_distribution'
include { PLOT_SCORE_DISTRIBUTION           } from '../modules/local/plot_score_distribution'
include { PREPARE_ENTITY_BINDING_RATIOS     } from '../modules/local/prepare_entity_binding_ratios'
include { PLOT_ENTITY_BINDING_RATIOS        } from '../modules/local/plot_entity_binding_ratios'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { PRODIGAL                    } from '../modules/nf-core/modules/prodigal/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METAPEP {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // // ####################################################################################################

    /*
    * Download proteins from entrez
    */
    DOWNLOAD_PROTEINS (
        INPUT_CHECK.out.ch_taxa_input
    )
    ch_versions = ch_versions.mix(DOWNLOAD_PROTEINS.out.versions)

    PRODIGAL(
        INPUT_CHECK.out.ch_nucl_input,
        "gff"
    )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    FRED2_GENERATEPEPTIDES (
        DOWNLOAD_PROTEINS.out.ch_entrez_fasta
        .mix (PRODIGAL.out.amino_acid_fasta)
    )
    ch_versions = ch_versions.mix(FRED2_GENERATEPEPTIDES.out.versions)

    SPLIT_PEPTIDES(
        FRED2_GENERATEPEPTIDES.out.splitted
    )
     ch_versions = ch_versions.mix(SPLIT_PEPTIDES.out.versions)

    // CREATE_PROTEIN_TSV (
    //     PRODIGAL.out.amino_acid_fasta
    // )
    // ch_versions = ch_versions.mix(CREATE_PROTEIN_TSV.out.versions)

    // /*
    //  * concat files and assign new, unique ids for all proteins (from different sources)
    //  */
    // GENERATE_PROTEIN_AND_ENTITY_IDS (
    //     CREATE_PROTEIN_TSV.out.ch_pred_proteins.collect { meta, file -> file }.ifEmpty([]),
    //     CREATE_PROTEIN_TSV.out.ch_pred_proteins.collect { meta, file -> meta }.ifEmpty([]),
    //     DOWNLOAD_PROTEINS.out.ch_entrez_proteins.ifEmpty([]),
    //     DOWNLOAD_PROTEINS.out.ch_entrez_entities_proteins.ifEmpty([]),
    //     DOWNLOAD_PROTEINS.out.ch_entrez_microbiomes_entities.ifEmpty([]),
    //     ch_proteins_input.collect { meta, file -> file }.ifEmpty([]),
    //     ch_proteins_input.collect { meta, file -> meta }.ifEmpty([])
    // )
    // ch_versions = ch_versions.mix(GENERATE_PROTEIN_AND_ENTITY_IDS.out.versions)

    // /*
    // * Generate peptides
    // */
    // GENERATE_PEPTIDES (
    //     GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_proteins
    // )
    // ch_versions = ch_versions.mix(GENERATE_PEPTIDES.out.versions)

    // /*
    // * Split prediction tasks (peptide, allele) into chunks of peptides that are to
    // * be predicted against the same allele for parallel prediction
    // */
    // SPLIT_PRED_TASKS (
    // GENERATE_PEPTIDES.out.ch_peptides,
    // GENERATE_PEPTIDES.out.ch_proteins_peptides,
    // GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities_proteins,
    // GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_microbiomes_entities_noweights,
    // INPUT_CHECK.out.ch_conditions,
    // INPUT_CHECK.out.ch_conditions_alleles,
    // INPUT_CHECK.out.ch_alleles
    // )
    // ch_versions = ch_versions.mix(SPLIT_PRED_TASKS.out.versions)

    // /*
    // * Perform epitope prediction
    // */
    // PREDICT_EPITOPES (
    //     SPLIT_PRED_TASKS.out.ch_epitope_prediction_chunks.flatten()
    // )
    // ch_versions = ch_versions.mix(PREDICT_EPITOPES.out.versions)

    // /*
    // * Merge prediction results from peptide chunks into one prediction result
    // */
    // // gather chunks of predictions and merge them already to avoid too many input files for `merge_predictions` process
    // // (causing "sbatch: error: Batch job submission failed: Pathname of a file, directory or other parameter too long")
    // // sort and buffer to ensure resume will work (inefficient, since this causes waiting for all predictions)
    // ch_epitope_predictions_buffered = PREDICT_EPITOPES.out.ch_epitope_predictions.toSortedList().flatten().buffer(size: 1000, remainder: true)
    // ch_epitope_prediction_warnings_buffered = PREDICT_EPITOPES.out.ch_epitope_prediction_warnings.toSortedList().flatten().buffer(size: 1000, remainder: true)

    // MERGE_PREDICTIONS_BUFFER (
    //     ch_epitope_predictions_buffered,
    //     ch_epitope_prediction_warnings_buffered
    // )
    // ch_versions = ch_versions.mix(MERGE_PREDICTIONS_BUFFER.out.versions)

    // MERGE_PREDICTIONS (
    //     MERGE_PREDICTIONS_BUFFER.out.ch_predictions_merged_buffer.collect(),
    //     MERGE_PREDICTIONS_BUFFER.out.ch_prediction_warnings_merged_buffer.collect()
    // )
    // ch_versions = ch_versions.mix(MERGE_PREDICTIONS.out.versions)

    // ASSIGN_NUCL_ENTITY_WEIGHTS (
    //     INPUT_CHECK.out.ch_weights,
    //     INPUT_CHECK.out.ch_conditions_weights
    // )
    // ch_versions = ch_versions.mix(ASSIGN_NUCL_ENTITY_WEIGHTS.out.versions)

    // /*
    //  * Create microbiome_entities
    //  */
    // FINALIZE_MICROBIOME_ENTITIES (
    //     DOWNLOAD_PROTEINS.out.ch_entrez_microbiomes_entities.ifEmpty([]),
    //     ASSIGN_NUCL_ENTITY_WEIGHTS.out.ch_nucl_microbiomes_entities.ifEmpty([]),
    //     GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_microbiomes_entities_noweights,
    //     GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities,
    //     INPUT_CHECK.out.ch_conditions
    // )
    // ch_versions = ch_versions.mix(FINALIZE_MICROBIOME_ENTITIES.out.versions)

    // /*
    // * Collect some numbers: proteins, peptides, unique peptides per conditon
    // */
    // if (!params.skip_stats){
    //     COLLECT_STATS (
    //         GENERATE_PEPTIDES.out.ch_peptides,
    //         GENERATE_PEPTIDES.out.ch_proteins_peptides,
    //         GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities_proteins,
    //         GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_microbiomes_entities_noweights,
    //         INPUT_CHECK.out.ch_conditions
    //     )
    //     ch_versions = ch_versions.mix(COLLECT_STATS.out.versions)
    // }

    // /*
    // * Generate figures
    // */
    // PREPARE_SCORE_DISTRIBUTION (
    //     MERGE_PREDICTIONS.out.ch_predictions,
    //     GENERATE_PEPTIDES.out.ch_proteins_peptides,
    //     GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities_proteins,
    //     FINALIZE_MICROBIOME_ENTITIES.out.ch_microbiomes_entities,
    //     INPUT_CHECK.out.ch_conditions,
    //     INPUT_CHECK.out.ch_conditions_alleles,
    //     INPUT_CHECK.out.ch_alleles
    // )
    // ch_versions = ch_versions.mix(PREPARE_SCORE_DISTRIBUTION.out.versions)

    // PLOT_SCORE_DISTRIBUTION (
    //     PREPARE_SCORE_DISTRIBUTION.out.ch_prep_prediction_scores.flatten(),
    //     INPUT_CHECK.out.ch_alleles,
    //     INPUT_CHECK.out.ch_conditions
    // )
    // ch_versions = ch_versions.mix(PLOT_SCORE_DISTRIBUTION.out.versions)

    // PREPARE_ENTITY_BINDING_RATIOS (
    //     MERGE_PREDICTIONS.out.ch_predictions,
    //     GENERATE_PEPTIDES.out.ch_proteins_peptides,
    //     GENERATE_PROTEIN_AND_ENTITY_IDS.out.ch_entities_proteins,
    //     FINALIZE_MICROBIOME_ENTITIES.out.ch_microbiomes_entities,
    //     INPUT_CHECK.out.ch_conditions,
    //     INPUT_CHECK.out.ch_conditions_alleles,
    //     INPUT_CHECK.out.ch_alleles
    // )
    // ch_versions = ch_versions.mix(PREPARE_ENTITY_BINDING_RATIOS.out.versions)

    // PLOT_ENTITY_BINDING_RATIOS (
    //     PREPARE_ENTITY_BINDING_RATIOS.out.ch_prep_entity_binding_ratios.flatten(),
    //     INPUT_CHECK.out.ch_alleles
    // )
    // ch_versions = ch_versions.mix(PLOT_ENTITY_BINDING_RATIOS.out.versions)

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    // //
    // // MODULE: MultiQC
    // //
    // workflow_summary    = WorkflowMetapep.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
    // ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
