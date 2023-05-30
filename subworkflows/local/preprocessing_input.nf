//
// Check input samplesheet and get read channels
//

include { INPUT_TO_DATAMODEL   } from '../../modules/local/input_to_datamodel'
include { UNPACK_BIN_ARCHIVES } from '../../modules/local/unpack_bin_archives'

workflow PREPROCESSING_INPUT {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    ch_versions = Channel.empty()

    INPUT_TO_DATAMODEL ( samplesheet )

    INPUT_TO_DATAMODEL.out.microbiomes
        // Read microbiomes table
        .splitCsv(sep:'\t', header:true)
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
        // TODO co-assembly case needs to be solved
        ch_microbiomes_branch.proteins
            .map { row ->
                    def meta = [:]
                    meta.id = row.microbiome_id
                    return [ meta, row.microbiome_path ]
                }
            .set { ch_proteins_input }
        ch_proteins_input.dump(tag:"proteins")

        // ASSEMBLY
            // Using the microbiome_bare_id to handle co-assembled input
            // microbiome_bare_id will be identical to microbiome_id if not co-assembled
            // The change in ID prevents redundant processes in protein prediction
        ch_microbiomes_branch.assembly
            .map { row ->
                    def meta = [:]
                    meta.id = row.microbiome_bare_id
                    meta.bin_basename = false
                    return [ meta, row.microbiome_path ]
                }
            .set { ch_assembly_input }
        ch_assembly_input.dump(tag:"assembly")

        // BINS
        // Using the microbiome_bare_id as meta.id to prevent redundant processing
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
                        meta.id = row.microbiome_bare_id
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
                    meta.id = row.microbiome_bare_id
                    return [ meta, row.microbiome_path ]
                }
            .set{ ch_microbiomes_bins_archives_packed }
            ch_microbiomes_bins_archives_packed.dump(tag:"bins_archives")

        /*
        * Unpack archived assembly bins
        */
        UNPACK_BIN_ARCHIVES(
            ch_microbiomes_bins_archives_packed.unique()
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

        // Concatenate the channels and remove redundant entries for nucleotide based inputs
        // In case of co-assembly the input fasta will be used for prediction only once
        ch_nucl_input           = ch_assembly_input.concat(ch_bins_archives_input, ch_bins_folders_input).unique()
        ch_nucl_input.dump(tag:"nucl")

        // ####################################################################################################

        ch_weights = Channel.empty()
        INPUT_TO_DATAMODEL.out.microbiomes
            .splitCsv(sep:'\t', header:true)
            .map { row ->
                    def meta = [:]
                    meta.id = row.microbiome_id
                    if (row.microbiome_type != 'taxa' && row.weights_path) [meta, row.weights_path]
                }
            .set { ch_weights }
            ch_weights.dump(tag:"weights")

    emit:
    ch_taxa_input
    ch_proteins_input
    ch_nucl_input
    ch_weights
    ch_microbiomes          = INPUT_TO_DATAMODEL.out.microbiomes
    ch_conditions           = INPUT_TO_DATAMODEL.out.conditions
    ch_alleles              = INPUT_TO_DATAMODEL.out.alleles
    ch_conditions_alleles   = INPUT_TO_DATAMODEL.out.conditions_alleles
    versions                = ch_versions
}
