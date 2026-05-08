/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HUMANN } from '../../modules/local/humann/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: FUNCTIONAL_ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This subworkflow performs functional annotation of metagenomic samples
    using HUMAnN3.

    DEPENDENCIES:
        - REQUIRES trimmed reads from fastp
        - REQUIRES MetaPhlAn4.1 output (taxonomic profile)
        - REQUIRES HUMAnN nucleotide database 
        - REQUIRES HUMAnN protein database

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_trimmed_reads        // channel: [meta, reads] - trimmed reads from fastp
    ch_metaphlan_profile    // channel: [meta, profile] - taxonomic profile from MetaPhlAn3
    ch_humann_nucleotide_db // path: path to HUMAnN nucleotide database
    ch_humann_protein_db    // path: path to HUMAnN protein database

    main:

    // Initialize ordered version tracking
    def ch_versions = channel.empty()

    // Prepare HUMAnN input by joining trimmed reads with MetaPhlAn profiles
    def ch_humann_input = ch_trimmed_reads
        .join(ch_metaphlan_profile)
        .map { meta, reads, profile -> [meta, reads, profile] }

    //
    // MODULE: Run HUMAnN for functional annotation
    //
    HUMANN(
        ch_humann_input,
        ch_humann_nucleotide_db,
        ch_humann_protein_db
    )
    ch_versions = ch_versions.concat(HUMANN.out.versions)

    emit:
    pathabundance   = HUMANN.out.pathabundance    // channel: [meta, pathabundance]
    pathcoverage    = HUMANN.out.pathcoverage     // channel: [meta, pathcoverage]
    genefamilies    = HUMANN.out.genefamilies     // channel: [meta, genefamilies] (if available)
    versions        = ch_versions                 // channel: versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
