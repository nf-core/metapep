/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetapep.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { UNPACK_BIN_ARCHIVES   } from '../modules/local/unpack_bin_archives'
include { DOWNLOAD_PROTEINS     } from '../modules/local/download_proteins'
include { CREATE_PROTEIN_TSV    } from '../modules/local/create_protein_tsv'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP                      } from '../modules/nf-core/modules/gunzip/main'
include { PRODIGAL                    } from '../modules/nf-core/modules/prodigal/main'
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
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

    INPUT_CHECK.out.ch_microbiomes
    // Read microbiomes table
    .splitCsv(header:true)
    // Convert paths to files
    .map {
        row ->
        row.microbiome_path = file(row.microbiome_path, checkIfExists: true)
        row
    }
    // Split into types
    .branch {
        row->
        taxa:      row.microbiome_type == 'taxa'
        proteins : row.microbiome_type == 'proteins'
        assembly:  row.microbiome_type == 'assembly'
        bins:      row.microbiome_type == 'bins'
    }
    .set{ ch_microbiomes_branch }

    // TAXA
    ch_microbiomes_branch.taxa
        .map { row ->
                def meta = [:]
                meta.id = row.microbiome_id
                return [ meta, row.microbiome_path ]
            }
        .set { ch_taxa_input }
    ch_taxa_input.dump(tag:"taxa")

    // PROTEINS
    ch_microbiomes_branch.proteins
        .map { row ->
                def meta = [:]
                meta.id = row.microbiome_id
                return [ meta, row.microbiome_path ]
            }
        .set { ch_proteins_input }
    ch_proteins_input.dump(tag:"proteins")

    // ASSEMBLY
    ch_microbiomes_branch.assembly
        .map { row ->
                def meta = [:]
                meta.id = row.microbiome_id
                meta.bin_basename = false
                return [ meta, row.microbiome_path ]
            }
        .set { ch_assembly_input }
    ch_assembly_input.dump(tag:"assembly")

    // BINS
    ch_microbiomes_branch.bins
        .branch {
                row ->
                folders : row.microbiome_path.isDirectory()
                archives : row.microbiome_path.isFile()
                other: true
            }
        .set{ ch_microbiomes_bins }

    // The file ending we expect for FASTA files
    fasta_suffix = ~/(?i)[.]fa(sta)?(.gz)?$/

    // BINS - LOCAL FOLDERS
    ch_microbiomes_bins.folders
        .flatMap { row ->
                def bin_files = row.microbiome_path.listFiles().findAll{ it.name =~ fasta_suffix }
                return bin_files.collect {
                    def meta = [:]
                    meta.id = row.microbiome_id
                    meta.bin_basename = it.name - fasta_suffix 
                    return [ meta, it ]
                }
            }
        .set { ch_bins_folders_input }
        ch_bins_folders_input.dump(tag:"bins_folders")

    // BINS - LOCAL OR REMOTE ARCHIVES
    ch_microbiomes_bins.archives
        .map { row ->
                def meta = [:]
                meta.id = row.microbiome_id
                return [ meta, row.microbiome_path ]
            }
        .set{ ch_microbiomes_bins_archives_packed }
        ch_microbiomes_bins_archives_packed.dump(tag:"bins_archives")
    
    /*
    * Unpack archived assembly bins
    */
    UNPACK_BIN_ARCHIVES(
        ch_microbiomes_bins_archives_packed
    )
    ch_versions = ch_versions.mix(UNPACK_BIN_ARCHIVES.out.versions)
    UNPACK_BIN_ARCHIVES.out.ch_microbiomes_bins_archives_unpacked.dump(tag:"bins_archives")

    ch_bins_archives_input = Channel.empty()
    UNPACK_BIN_ARCHIVES.out.ch_microbiomes_bins_archives_unpacked
        .flatMap { meta, bin_files ->
                bin_files = bin_files.findAll{ it.name =~ fasta_suffix }
                if (bin_files.isEmpty()) log.warn("WARNING - Archive provided for microbiome ID ${microbiome_id} did not yield any bin files")
                return bin_files.collect {
                    def meta_new = [:]
                    meta_new.id = meta.id
                    meta_new.bin_basename = it.name - fasta_suffix
                    return [ meta_new, it ]
                }
            }
        .set{ ch_bins_archives_input }
    ch_bins_archives_input.dump(tag:"bins_archives")

    // Concatenate the channels for nucleotide based inputs
    ch_nucl_input           = ch_assembly_input.concat(ch_bins_archives_input, ch_bins_folders_input)
    ch_nucl_input.dump(tag:"nucl")

    // ####################################################################################################

    ch_weights = Channel.empty()
    INPUT_CHECK.out.ch_microbiomes
        .splitCsv(header:true)
        .map { row ->
                def meta = [:]
                meta.id = row.microbiome_id
                if (row.microbiome_type != 'taxa' && row.weights_path) [meta, row.weights_path]
            }
        .set { ch_weights }
        ch_weights.dump(tag:"weights")

    /*
    * Download proteins from entrez
    */
    DOWNLOAD_PROTEINS (
        ch_taxa_input.map { meta, file -> meta.id }.collect().dump(tag:"taxa"),
        ch_taxa_input.map { meta, file -> file }.collect().dump(tag:"taxa")
    )
    ch_versions = ch_versions.mix(DOWNLOAD_PROTEINS.out.versions)

    ch_nucl_input
        .branch{ meta, file -> 
            zipped: file.name =~ ~/(?i)[.]gz$/ 
            unzipped: true
        }
        .set {ch_nucl_unzip}

    GUNZIP(
        ch_nucl_unzip.zipped
    )
    ch_versions = ch_versions.mix(GUNZIP.out.versions)

    ch_nucl_input_unzipped = GUNZIP.out.gunzip.concat(ch_nucl_unzip.unzipped)

    PRODIGAL(
        ch_nucl_input_unzipped,
        "gff"
    )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    CREATE_PROTEIN_TSV (
        PRODIGAL.out.amino_acid_fasta
    )
    ch_versions = ch_versions.mix(CREATE_PROTEIN_TSV.out.versions)

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
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
