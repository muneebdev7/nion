/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METAPHLAN3_METAPHLAN3  } from '../../modules/nf-core/metaphlan3/metaphlan3/main'
include { METAPHLAN3_MERGEMETAPHLANTABLES } from '../../modules/nf-core/metaphlan3/mergemetaphlantables/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: TAXONOMIC_CLASSIFICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This subworkflow performs taxonomic classification and abundance estimation
    using MetaPhlAn.

    DEPENDENCIES:
        - REQUIRES trimmed reads from fastp
        - REQUIRES MetaPhlAn database

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TAXONOMIC_CLASSIFICATION {

    take:
    ch_trimmed_reads    // channel: [meta, reads] - trimmed reads from fastp
    ch_metaphlan_db     // path: path to MetaPhlAn3 database

    main:

    // Initialize ordered version tracking
    def ch_versions = channel.empty()

    //
    // MODULE: Run MetaPhlAn3 for taxonomic classification
    //
    METAPHLAN3_METAPHLAN3(
        ch_trimmed_reads.map { meta, reads -> [meta, reads] },
        ch_metaphlan_db
    )
    ch_versions = ch_versions.concat(METAPHLAN3_METAPHLAN3.out.versions)

    //
    // MODULE: Merge MetaPhlAn profiles across all samples
    //
    METAPHLAN3_MERGEMETAPHLANTABLES(
        METAPHLAN3_METAPHLAN3.out.profile
            .collect { _meta, profile -> profile }
            .map { profiles -> [[id:'merged'], profiles] }
    )
    ch_versions = ch_versions.concat(METAPHLAN3_MERGEMETAPHLANTABLES.out.versions)

    emit:
    profile         = METAPHLAN3_METAPHLAN3.out.profile             // channel: [meta, profile]
    merged_table    = METAPHLAN3_MERGEMETAPHLANTABLES.out.txt       // channel: [meta, merged_table]
    versions        = ch_versions                                   // channel: versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
