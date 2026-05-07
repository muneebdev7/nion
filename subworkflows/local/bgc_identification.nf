/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP as GUNZIP_CONTIGS  } from '../../modules/nf-core/gunzip/main'
include { ANTISMASH_ANTISMASH       } from '../../modules/nf-core/antismash/antismash/main'
include { COVERM_CONTIG             } from '../../modules/nf-core/coverm/contig/main'
include { COVERM_MERGE_TABLES       } from '../../modules/local/coverm_merge_tables/main'
include { ANTISMASH_PROCESS_SAMPLES } from '../../modules/local/antismash_process_samples/main'
include { ANTISMASH_EXTRACT_BGCS    } from '../../modules/local/antismash_extract_bgcs/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BGC_IDENTIFICATION 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This subworkflow performs biosynthetic gene cluster (BGC) identification
    using antiSMASH on assembled contigs.

    DEPENDENCIES:
        - REQUIRES trimmed reads from fastp
        - REQUIRES assembled contigs from MEGAHIT

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BGC_IDENTIFICATION {

    take:
    ch_reads   // channel: [meta, reads] - trimmed reads from fastp
    ch_contigs // channel: [meta, contigs] - assembled contigs from MEGAHIT

    main:

    // Initialize ordered version tracking
    def ch_versions = channel.empty()
    def ch_antismash_databases = channel.fromPath(
        params.antismash_databases,
        checkIfExists: true
    ).first()
    def gff_input = params.gff ? file(params.gff, checkIfExists: true) : ''

    // MODULE 1: Unzip MEGAHIT contigs
    GUNZIP_CONTIGS(ch_contigs)
    def ch_contigs_unzipped = GUNZIP_CONTIGS.out.gunzip

    // MODULE 2: Run antiSMASH for BGC identification
    ANTISMASH_ANTISMASH(
        ch_contigs_unzipped,
        ch_antismash_databases,
        channel.value(gff_input)
    )
    ch_versions = ch_versions.concat(ANTISMASH_ANTISMASH.out.versions)

    // MODULE 3: Process antiSMASH GBK outputs to build regions TSVs
    ANTISMASH_PROCESS_SAMPLES(
        ANTISMASH_ANTISMASH.out.gbk_input,
    )

    // MODULE 4: Extract BGC FASTAs from contigs using regions TSVs
    ANTISMASH_EXTRACT_BGCS(
        ch_contigs_unzipped.join(ANTISMASH_PROCESS_SAMPLES.out.regions, by: 0)
            .map { meta, contigs, regions_tsv -> [meta, contigs, regions_tsv] },
    )

    // MODULE 5: Run CoverM contig coverage after antiSMASH
    def ch_coverm_input = ch_reads
        .join(ANTISMASH_EXTRACT_BGCS.out.bgc_fasta, by: 0)
        .map { meta, reads, bgc_fasta -> [meta, reads, bgc_fasta] }

    COVERM_CONTIG(
        ch_coverm_input.map { meta, reads, _bgc_fasta -> [meta, reads] },
        ch_coverm_input.map { meta, _reads, bgc_fasta -> [meta, bgc_fasta] },
        channel.value(false),
        channel.value(false),
        channel.value(false)
    )

    // MODULE 6: Merge CoverM depth tables and standardize headers
    COVERM_MERGE_TABLES(
        COVERM_CONTIG.out.coverage
            .collect { _meta, coverage -> coverage }
            .map { coverage_files -> [[id: 'coverm_merged'], coverage_files] }
    )

    emit:
    html            = ANTISMASH_ANTISMASH.out.html          // channel: [meta, index.html]
    gbk_results     = ANTISMASH_ANTISMASH.out.gbk_results   // channel: [meta, region gbk]
    json_results    = ANTISMASH_ANTISMASH.out.json_results  // channel: [meta, json]
    regions_tsv     = ANTISMASH_PROCESS_SAMPLES.out.regions // channel: [meta, regions tsv]
    bgc_fasta       = ANTISMASH_EXTRACT_BGCS.out.bgc_fasta  // channel: [meta, bgc fasta]
    coverm_coverage = COVERM_CONTIG.out.coverage            // channel: [meta, *.depth.tsv]
    coverm_bam      = COVERM_CONTIG.out.bam                 // channel: [meta, *.bam]
    coverm_merged   = COVERM_MERGE_TABLES.out.merged        // channel: [meta, merged coverm table]
    versions        = ch_versions                           // channel: versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
