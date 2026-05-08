/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METAPHLAN4              } from '../../modules/local/metaphlan4/main'
include { METAPHLAN4_MERGE_TABLES } from '../../modules/local/metaphlan4_merge_tables/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: TAXONOMIC_CLASSIFICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This subworkflow performs taxonomic classification and abundance estimation
    using MetaPhlAn4.

    DEPENDENCIES:
        - REQUIRES trimmed reads from fastp
        - REQUIRES MetaPhlAn4.1 database index named as 'mpa_vJun23_CHOCOPhlAnSGB_202307'

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TAXONOMIC_CLASSIFICATION {

    take:
    ch_trimmed_reads    // channel: [meta, reads] - trimmed reads from fastp
    ch_metaphlan_db     // path: path to MetaPhlAn4.1 database

    main:

    // Initialize ordered version tracking
    def ch_versions = channel.empty()

    //
    // MODULE: Run MetaPhlAn4 for taxonomic classification
    //
    METAPHLAN4(
        ch_trimmed_reads.map { meta, reads -> [meta, reads] },
        ch_metaphlan_db
    )
    ch_versions = ch_versions.concat(METAPHLAN4.out.versions)

    //
    // MODULE: Merge MetaPhlAn profiles across all samples
    //
    METAPHLAN4_MERGE_TABLES(
        METAPHLAN4.out.profile
            .collect { _meta, profile -> profile }
            .map { profiles -> [[id:'merged'], profiles] }
    )
    ch_versions = ch_versions.concat(METAPHLAN4_MERGE_TABLES.out.versions)

    emit:
    profile         = METAPHLAN4.out.profile            // channel: [meta, profile]
    merged_table    = METAPHLAN4_MERGE_TABLES.out.txt   // channel: [meta, merged_table]
    versions        = ch_versions                       // channel: versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
