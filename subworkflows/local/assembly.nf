/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MEGAHIT } from '../../modules/nf-core/megahit/main'



    /*
    ================================================================================
                                    Assembly
    ================================================================================
    */

workflow ASSEMBLY {

    take:
    ch_trimmed_reads    // channel: [meta, reads] - trimmed reads from fast

    main:
    // Initialize ordered version tracking
    def ch_versions = channel.empty()

    // Prepare FASTP output for MEGAHIT input
    def ch_megahit_input = ch_trimmed_reads.map { meta, reads ->
        def reads1 = meta.single_end ? reads : reads[0]
        def reads2 = meta.single_end ? [] : reads[1]
        [meta, reads1, reads2]
    }

    //
    // MODULE: Run MEGAHIT for assembly
    //
    MEGAHIT (
        ch_megahit_input
    )
    // Collect version information
    ch_versions = ch_versions.concat(MEGAHIT.out.versions)

    emit:
        contigs         = MEGAHIT.out.contigs   // channel: [meta, contigs]
        versions        = ch_versions           // channel: versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
