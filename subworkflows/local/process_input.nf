//
// Check input samplesheet, get main data channels and create data tables with provided data (i.e. microbiomes, conditions, alleles, conditions_alleles)
//

include { CHECK_SAMPLESHEET_CREATE_TABLES   } from '../../modules/local/check_samplesheet_create_tables'
include { UNPACK_BIN_ARCHIVES               } from '../../modules/local/unpack_bin_archives'
include { UNIFY_MODEL_LENGTHS               } from '../../modules/local/unify_model_lengths'

workflow PROCESS_INPUT {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    ch_versions = Channel.empty()

    CHECK_SAMPLESHEET_CREATE_TABLES ( samplesheet )
    ch_versions = ch_versions.mix(CHECK_SAMPLESHEET_CREATE_TABLES.out.versions)

    // When PSSMs methods are used we can check in epytope if the model exists for the given peptide lengths
    // Only intersection of allele model lengths are used in further analysis
    if (params.pred_method == "syfpeithi" && !params.allow_inconsistent_pep_lengths) {
        UNIFY_MODEL_LENGTHS (CHECK_SAMPLESHEET_CREATE_TABLES.out.samplesheet_valid)
        peptide_lengths =  UNIFY_MODEL_LENGTHS.out.unified_pep_lens.first().tokenize(",")

        // Throw Error or warning depending on file prefix
        log_prefix = UNIFY_MODEL_LENGTHS.out.log.first().baseName.tokenize("_").flatten().first()
        log_prefix.map{it ->
            if (it == "WARNING"){
                log.warn "There is no SYFPEITHI model available for at least one allele and the peptide lengths are reduced to the common denominator of the alleles."
                log.warn "For more information about which peptide lengths were used check '$params.outdir/logs/WARNING_unify_peptide_lengths.log'"
            } else if (it == "ERROR") {
                error "No models with matching peptide lengths could be found. Check '$params.outdir/logs/ERROR_unify_peptide_lengths.log' for further details."
            }
        }
        ch_versions = ch_versions.mix(UNIFY_MODEL_LENGTHS.out.versions)
    } else {
        // for MHCFlurry and MHCnuggets assing the input lengths
        peptide_lengths = Channel.fromList( params.min_pep_len..params.max_pep_len )
    }

    CHECK_SAMPLESHEET_CREATE_TABLES.out.microbiomes
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
        CHECK_SAMPLESHEET_CREATE_TABLES.out.microbiomes
            .splitCsv(sep:'\t', header:true)
            .map { row ->
                    def meta = [:]
                    meta.id = row.microbiome_id
                    if (row.microbiome_type != 'taxa' && row.weights_path) [meta, row.weights_path]
                }
            .set { ch_weights }
            ch_weights.dump(tag:"weights")

    CHECK_SAMPLESHEET_CREATE_TABLES.out.alleles
        .countLines()
        .map{
            it ->
                if (it.toInteger()-1 > params.max_task_num) {
                    error "Number of chosen alleles is larger than the number of maximum tasks submitted by the peptide epitope prediction process. For each allele at least 1 task needs to be executed."
                }
            }

    emit:
    peptide_lengths         = peptide_lengths.collect()
    ch_taxa_input
    ch_nucl_input
    ch_weights
    ch_microbiomes          = CHECK_SAMPLESHEET_CREATE_TABLES.out.microbiomes
    ch_conditions           = CHECK_SAMPLESHEET_CREATE_TABLES.out.conditions
    ch_alleles              = CHECK_SAMPLESHEET_CREATE_TABLES.out.alleles
    ch_conditions_alleles   = CHECK_SAMPLESHEET_CREATE_TABLES.out.conditions_alleles
    versions                = ch_versions
}
