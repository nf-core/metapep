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
include { FRED2_GENERATEPEPTIDES            } from '../modules/local/fred2_generatepeptides.nf'
include { SORT_PEPTIDES                     } from '../modules/local/sort_peptides'
include { REMOVE_DUPLICATE_PEPTIDES         } from '../modules/local/remove_duplicate_peptides'
include { SPLIT_PEPTIDES                    } from '../modules/local/split_peptides.nf'
include { GET_PREDICTION_VERSIONS           } from '../modules/local/get_prediction_versions'
include { PEPTIDE_PREDICTION                } from '../modules/local/peptide_prediction'
include { CAT_FILES as CAT_TSV              } from '../modules/local/cat_files'
include { CAT_FILES as CAT_FASTA            } from '../modules/local/cat_files'
include { CSVTK_CONCAT                      } from '../modules/local/csvtk_concat'
include { MERGE_JSON as MERGE_JSON_SINGLE   } from '../modules/local/merge_json'
include { MERGE_JSON as MERGE_JSON_MULTI    } from '../modules/local/merge_json'
include { PREPARE_PLOTS                     } from '../modules/local/prepare_plots'
include { PLOT_SCORE_DISTRIBUTION           } from '../modules/local/plot_score_distribution'
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

    // ####################################################################################################
    /*
    * Download proteins from entrez
    */
    DOWNLOAD_PROTEINS (
        INPUT_CHECK.out.ch_taxa_input
    )
    ch_versions = ch_versions.mix(DOWNLOAD_PROTEINS.out.versions)

    /*
    * Predict proteins using prodigal
    */
    PRODIGAL(
        INPUT_CHECK.out.ch_nucl_input,
        "gff"
    )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    /*
    * Generate peptides from proteins
    */
    FRED2_GENERATEPEPTIDES (
        DOWNLOAD_PROTEINS.out.ch_entrez_fasta
        .mix (PRODIGAL.out.amino_acid_fasta)
        .splitFasta( by: params.max_fasta_size, file:true )
    )
    ch_versions = ch_versions.mix(FRED2_GENERATEPEPTIDES.out.versions)

    /*
    * Sort chunks of peptides alphabetically to enable external merge sort of process REMOVE_DUPLICATE_PEPTIDES
    */
    SORT_PEPTIDES (
        FRED2_GENERATEPEPTIDES.out.splitted
    )

    SORT_PEPTIDES.out.output
    .map { meta, file ->
        [meta.alleles.tokenize(';'), meta, file]
    }
    .transpose()
    .groupTuple()
    .set {ch_sorted_peptides}

    /*
    * External merge sorting of peptides, removing peptides in the process
    * Writing additional info to csv that is given to PEPTIDE_PREDICTION
    */
    REMOVE_DUPLICATE_PEPTIDES (
        ch_sorted_peptides,
        INPUT_CHECK.out.ch_weights.dump(tag:'weights')
    )

    /*
    * Split peptides into chunks for parallel prediction
    */
    SPLIT_PEPTIDES(
        REMOVE_DUPLICATE_PEPTIDES.out.output.map { allele, file ->
        meta = [:]
        meta.alleles = allele
        meta.sample = allele
        return [meta, file]
    }
    )
    ch_versions = ch_versions.mix(SPLIT_PEPTIDES.out.versions)

    /*
    * Acquire versions of prediction tools
    */
    GET_PREDICTION_VERSIONS([])
    ch_prediction_tool_versions = GET_PREDICTION_VERSIONS.out.versions.ifEmpty("")

    /*
    * Predict epitopes
    */
    PEPTIDE_PREDICTION (
        SPLIT_PEPTIDES
            .out
            .splitted
            .combine( ch_prediction_tool_versions )
            .transpose()
    )
    ch_versions = ch_versions.mix( PEPTIDE_PREDICTION.out.versions )

    PEPTIDE_PREDICTION
        .out
        .predicted
        .groupTuple()
        .map { meta, predicted ->
            meta.files = predicted.size()
            return [meta, predicted]}
        .branch {
            meta_data, predicted ->
                multi: meta_data.files > 1
                    return [ meta_data, predicted ]
                single: meta_data.files == 1
                    return [ meta_data, predicted ]
        }
        .set { ch_predicted_peptides }

    /*
    * Combine epitope prediction results
    */
    CAT_TSV(
        ch_predicted_peptides.single
    )
    CSVTK_CONCAT(
        ch_predicted_peptides.multi
    )
    ch_versions = ch_versions.mix( CSVTK_CONCAT.out.versions)

    ch_predictions = CAT_TSV.out.output.mix( CSVTK_CONCAT.out.predicted )

    /*
    * Combine protein sequences (if they exist)
    */
    CAT_FASTA(
        PEPTIDE_PREDICTION
            .out
            .fasta
            .groupTuple()
    )

    PEPTIDE_PREDICTION
        .out
        .json
        .groupTuple()
        .map { meta, json ->
            meta.files = json.size()
            return [meta, json]}
        .branch {
            meta, json ->
                multi: meta.files > 1
                    return [ meta, json ]
                single: meta.files == 1
                    return [ meta, json ]
        }
        .set { ch_json_reports }

    /*
    * Combine epitope prediction reports
    */
    MERGE_JSON_SINGLE(
        ch_json_reports.single
    )
    MERGE_JSON_MULTI(
        ch_json_reports.multi
    )
    ch_versions = ch_versions.mix( MERGE_JSON_SINGLE.out.versions.ifEmpty(null) )
    ch_versions = ch_versions.mix( MERGE_JSON_MULTI.out.versions.ifEmpty(null) )

    /*
    * Create tables used for plotting results
    */
    PREPARE_PLOTS(
        ch_predictions
    )
    ch_versions = ch_versions.mix(PREPARE_PLOTS.out.versions)

    /*
    * Plot distribution of prediction scores
    */
    PLOT_SCORE_DISTRIBUTION (
        PREPARE_PLOTS.out.ch_prep_prediction_scores
    )
    ch_versions = ch_versions.mix(PLOT_SCORE_DISTRIBUTION.out.versions)

    /*
    * Plot entity binding ratios
    */
    PLOT_ENTITY_BINDING_RATIOS (
        PREPARE_PLOTS.out.ch_prep_entity_binding_ratios
    )
    ch_versions = ch_versions.mix(PLOT_ENTITY_BINDING_RATIOS.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMetapep.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
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
