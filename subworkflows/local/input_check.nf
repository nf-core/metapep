//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )

    emit:
    ch_microbiomes          = SAMPLESHEET_CHECK.out.microbiomes
    ch_conditions           = SAMPLESHEET_CHECK.out.conditions
    ch_alleles              = SAMPLESHEET_CHECK.out.alleles
    ch_conditions_alleles   = SAMPLESHEET_CHECK.out.conditions_alleles
    versions                = SAMPLESHEET_CHECK.out.versions            // channel: [ versions.yml ]
}
