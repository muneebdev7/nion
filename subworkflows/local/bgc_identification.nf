/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANTISMASH_ANTISMASH } from '../../modules/nf-core/antismash/antismash/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BGC_IDENTIFICATION 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This subworkflow performs biosynthetic gene cluster (BGC) identification
    using antiSMASH on assembled contigs.

    DEPENDENCIES:
        - REQUIRES trimmed reads from fastp
        - REQUIRES MetaPhlAn3 output (taxonomic profile)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BGC_IDENTIFICATION {

    take:
    ch_contigs // channel: [meta, contigs] - assembled contigs from MEGAHIT

    main:

    // Initialize ordered version tracking
    def ch_versions = channel.empty()
    def ch_antismash_databases = channel.fromPath(
        params.antismash_databases,
        checkIfExists: true
    ).first()
    def gff_input = params.gff ? file(params.gff, checkIfExists: true) : ''

    // MODULE: Run antiSMASH for BGC identification
    ANTISMASH_ANTISMASH(
        ch_contigs,
        ch_antismash_databases,
        channel.value(gff_input)
    )
    ch_versions = ch_versions.concat(ANTISMASH_ANTISMASH.out.versions)

    emit:
    html            = ANTISMASH_ANTISMASH.out.html       // channel: [meta, index.html]
    gbk_results     = ANTISMASH_ANTISMASH.out.gbk_results // channel: [meta, region gbk]
    json_results    = ANTISMASH_ANTISMASH.out.json_results // channel: [meta, json]
    versions        = ch_versions                         // channel: versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
