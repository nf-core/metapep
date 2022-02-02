//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK     } from '../../modules/local/samplesheet_check'
include { UNPACK_BIN_ARCHIVES   } from '../../modules/local/unpack_bin_archives'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    ch_versions = Channel.empty()

    SAMPLESHEET_CHECK ( samplesheet )
        .microbiomes
        .splitCsv ( header:true, sep:'\t' )
        .map { create_meta(it) }
        .branch {
        meta, file->
        taxa:           meta.type == 'taxa'
        proteins :      meta.type == 'proteins'
        assembly:       meta.type == 'assembly'
        bins_folders :  meta.type == 'bins' && file.isDirectory()
        bins_archives : meta.type == 'bins' && file.isFile()
        }
        .set{ ch_microbiomes_branch }
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    /*
    * Unpack archived assembly bins
    */
    UNPACK_BIN_ARCHIVES(
        ch_microbiomes_branch.bins_archives
    )
    ch_versions = ch_versions.mix(UNPACK_BIN_ARCHIVES.out.versions)

    // The file ending we expect for FASTA files
    fasta_suffix = ~/(?i)[.]fa(sta)?(.gz)?$/

    UNPACK_BIN_ARCHIVES.out.ch_microbiomes_bins_archives_unpacked.concat(ch_microbiomes_branch.bins_folders).dump(tag:"test")
        .flatMap { meta, bin_files ->
                bin_files = bin_files instanceof List ? bin_files.findAll{ it.name =~ fasta_suffix } : bin_files.listFiles().findAll{ it.name =~ fasta_suffix }
                if (bin_files.isEmpty()) log.warn("WARNING - Archive or folder provided for microbiome ID ${meta.id} did not yield any bin files")
                return bin_files.collect {
                    def meta_new = [:]
                    meta_new.id             = meta.id
                    meta_new.sample         = meta.sample
                    meta_new.conditions     = meta.conditions
                    meta_new.alleles        = meta.alleles
                    meta_new.type           = meta.type
                    meta_new.bin_basename   = it.name - fasta_suffix
                    return [ meta_new, it ]
                }
            }
        .set{ ch_bins_input }

    // Concatenate the channels for nucleotide based inputs
    ch_nucl_input           = ch_microbiomes_branch.assembly.concat(ch_bins_input)
    ch_nucl_input.dump(tag:"nucl")

    ch_microbiomes_branch.taxa
        .map { meta, file -> [meta, file.splitCsv(sep:"\t", header:true)] }
        .transpose()
        .map { meta, taxon ->
            return [meta, taxon.taxon_id]}
        .groupTuple(by:1)
        .map { meta, taxon_id ->
            def meta_new = [:]
            meta_new.id             = taxon_id
            meta_new.sample         = taxon_id
            meta_new.alleles        = meta.collect {
                microbiome ->
                microbiome.alleles.split('[; ]')
            }
            .flatten()
            .unique().join(';')
            meta_new.microbiomes    = meta.collect { m ->
                def meta_reduced = [:]
                meta_reduced.id             = m.id
                meta_reduced.conditions     = m.conditions
                meta_reduced.type           = m.type
                meta_reduced.bin_basename   = m.bin_basename
                return meta_reduced
            }
            return [meta_new, taxon_id]
        }
        .set{ ch_taxa_input }

    emit:
    ch_taxa_input
    ch_proteins_input           = ch_microbiomes_branch.proteins
    ch_nucl_input
    ch_microbiomes              = SAMPLESHEET_CHECK.out.microbiomes
    ch_conditions               = SAMPLESHEET_CHECK.out.conditions
    ch_conditions_microbiomes   = SAMPLESHEET_CHECK.out.conditions_microbiomes
    ch_alleles                  = SAMPLESHEET_CHECK.out.alleles
    ch_conditions_alleles       = SAMPLESHEET_CHECK.out.conditions_alleles
    ch_weights                  = SAMPLESHEET_CHECK.out.weights
    ch_conditions_weights       = SAMPLESHEET_CHECK.out.conditions_weights
    versions                    = ch_versions                               // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_meta(LinkedHashMap row) {
    def meta = [:]
    meta.id             = row.microbiome_id
    meta.sample         = row.microbiome_id
    meta.conditions     = row.conditions
    meta.alleles        = row.alleles
    meta.type           = row.microbiome_type
    meta.bin_basename   = false
    return [ meta, file(row.microbiome_path, checkIfExists: true) ]
}
