/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_INDEX                            } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM                              } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX                       } from '../../modules/nf-core/samtools/index/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This subworkflow performs indexing, mapping (BWA index & BWA mem) and 
    sorting of alignments of metagenomic samples.

    DEPENDENCIES:
        - REQUIRES trimmed reads from fastp
        - REQUIRES assembled contigs from MEGAHIT

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAPPING {
    
    take:
    ch_trimmed_reads  // channel: [meta, reads] - cleaned/trimmed reads
    ch_assemblies     // channel: [meta, contigs] - assembled contigs

    main:

    // Initialize ordered version tracking
    def ch_versions = channel.empty()

    //
    // MODULE: Build BWA index for assembled contigs
    //
    BWA_INDEX(ch_assemblies)
    def ch_index_for_bwamem = BWA_INDEX.out.index.map { meta, index ->
        [meta, index]
    }

    //
    // MODULE: Run BWA_MEM for alignment of trimmed reads to indexed contigs
    //
    def ch_bwamem_input = ch_trimmed_reads
        .join(ch_index_for_bwamem, by: 0)
        .join(ch_assemblies, by: 0)
        .map { meta, reads, index, contigs ->
            [meta, reads, index, contigs]
        }
    
    BWA_MEM(
        ch_bwamem_input.map { meta, reads, _index, _contigs -> [meta, reads] },
        ch_bwamem_input.map { meta, _reads, index, _contigs -> [meta, index] },
        ch_bwamem_input.map { meta, _reads, _index, contigs -> [meta, contigs] },
        true
    )

    //
    // MODULE: Run SAMtools/Index to index the aligned BAM files
    //
    SAMTOOLS_INDEX(BWA_MEM.out.bam)

    // Output channels for downstream use
    def ch_sorted_bam = BWA_MEM.out.bam

    emit:
    sorted_bam = ch_sorted_bam
    versions   = ch_versions

}