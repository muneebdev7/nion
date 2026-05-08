/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP as GUNZIP_BINS                } from '../../modules/nf-core/gunzip/main'
include { SEQFU_DEREP                          } from '../../modules/nf-core/seqfu/derep/main'
include { GUNZIP as GUNZIP_DEREP               } from '../../modules/nf-core/gunzip/main'
include { GTDBTK_CLASSIFYWF                    } from '../../modules/nf-core/gtdbtk/classifywf/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BIN_TAXONOMIC_CLASSIFICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This subworkflow performs post-binning taxonomic classification.
    It decompresses MetaBAT2 bins, dereplicates them, then runs GTDB-Tk.

    DEPENDENCIES:
        - REQUIRES gzipped bins from BINNING (MetaBAT2)
        - REQUIRES GTDB-Tk database directory

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BIN_TAXONOMIC_CLASSIFICATION {

    take:
    ch_bins_gz       // channel: [meta, bin_fasta_gz]
    ch_gtdbtk_db     // channel: path to GTDB-Tk database root directory

    main:

    def ch_versions = channel.empty()

    // Decompress each gzipped bin produced by MetaBAT2
    GUNZIP_BINS(ch_bins_gz)

    // Group decompressed bins per sample and dereplicate within each sample
    def ch_bins_for_derep = GUNZIP_BINS.out.gunzip
        .groupTuple(by: 0)
        .map { meta, bins -> [meta, bins] }

    SEQFU_DEREP(ch_bins_for_derep)

    // GTDB-Tk classify_wf expects uncompressed genomes in genome_dir
    GUNZIP_DEREP(SEQFU_DEREP.out.fasta)

    def ch_gtdbtk_input = GUNZIP_DEREP.out.gunzip
        .map { meta, derep_fasta -> [meta, [derep_fasta]] }

    def ch_gtdbtk_db_input = ch_gtdbtk_db
        .map { db -> ['gtdbtk_db', db] }

    GTDBTK_CLASSIFYWF(
        ch_gtdbtk_input,
        ch_gtdbtk_db_input,
        params.gtdbtk_use_pplacer_scratch_dir
    )

    emit:
    bins_unzipped    = GUNZIP_BINS.out.gunzip
    derep_bins       = SEQFU_DEREP.out.fasta
    derep_unzipped   = GUNZIP_DEREP.out.gunzip
    gtdb_summary     = GTDBTK_CLASSIFYWF.out.summary
    gtdb_outdir      = GTDBTK_CLASSIFYWF.out.gtdb_outdir
    gtdb_tree        = GTDBTK_CLASSIFYWF.out.tree
    gtdb_markers     = GTDBTK_CLASSIFYWF.out.markers
    gtdb_msa         = GTDBTK_CLASSIFYWF.out.msa
    gtdb_user_msa    = GTDBTK_CLASSIFYWF.out.user_msa
    gtdb_filtered    = GTDBTK_CLASSIFYWF.out.filtered
    gtdb_failed      = GTDBTK_CLASSIFYWF.out.failed
    gtdb_log         = GTDBTK_CLASSIFYWF.out.log
    gtdb_warnings    = GTDBTK_CLASSIFYWF.out.warnings
    versions         = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/