/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_nion_pipeline'

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                    } from '../modules/nf-core/fastqc/main'
include { FASTP                     } from '../modules/nf-core/fastp/main'

//
// SUBWORKFLOWS
//
include { TAXONOMIC_CLASSIFICATION  } from '../subworkflows/local/taxonomic_classification'
include { FUNCTIONAL_ANNOTATION     } from '../subworkflows/local/functional_annotation'
include { ASSEMBLY                  } from '../subworkflows/local/assembly'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NION {

    take:
    ch_raw_short_reads // channel: samplesheet read in from --input
    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_raw_short_reads
    )
    ch_versions = ch_versions.concat(FASTQC.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, zip -> zip })

    //
    // MODULE: Run Fastp for trimming and filtering
    //
    FASTP (
        ch_raw_short_reads.map { meta, reads -> [meta, reads, []] },
        false,
        params.save_trimmed_fail,
        false
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { _meta, json -> json })

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Subworkflow 1: TAXONOMIC CLASSIFICATION - Taxonomic profiling using MetaPhlAn3
        OPTIONAL: Only runs if --metaphlan_db is provided
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.metaphlan_db) {
        def ch_metaphlan_db = channel.fromPath(
            params.metaphlan_db,
            checkIfExists: true
        ).first()

        TAXONOMIC_CLASSIFICATION(
            FASTP.out.reads,
            ch_metaphlan_db
        )

        ch_versions = ch_versions.concat(TAXONOMIC_CLASSIFICATION.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(TAXONOMIC_CLASSIFICATION.out.profile.map { _meta, profile -> profile })
    } else {
        log.warn("MetaPhlAn is disabled: provide --metaphlan_db.")
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Subworkflow 2: FUNCTIONAL ANNOTATION - Functional profiling using HUMAnN
        OPTIONAL: Only runs if both --humann_nucleotide_db and --humann_protein_db provided
        REQUIRED: --metaphlan_db must be enabled (MetaPhlAn output is needed as input)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.humann_nucleotide_db && params.humann_protein_db) {

        // Validate that MetaPhlAn is enabled (required for HUMAnN)
        if (!params.metaphlan_db) {
            log.error(
                "ERROR: --metaphlan_db is required when running HUMAnN.\n" +
                "       HUMAnN requires MetaPhlAn taxonomic output as input.\n" +
                "       Please provide --metaphlan_db or disable HUMAnN (omit --humann_nucleotide_db and --humann_protein_db)."
            )
            System.exit(1)
        }

        def ch_humann_nucleotide_db = channel.fromPath(
            params.humann_nucleotide_db,
            checkIfExists: true
        ).first()
        def ch_humann_protein_db = channel.fromPath(
            params.humann_protein_db,
            checkIfExists: true
        ).first()

        FUNCTIONAL_ANNOTATION(
            FASTP.out.reads,
            TAXONOMIC_CLASSIFICATION.out.profile,
            ch_humann_nucleotide_db,
            ch_humann_protein_db
        )

        ch_versions = ch_versions.concat(FUNCTIONAL_ANNOTATION.out.versions)
    } else {
        if (!params.humann_nucleotide_db || !params.humann_protein_db) {
            log.info(
                "HUMAnN is disabled. To enable, provide both:\n" +
                "  --humann_nucleotide_db <path>\n" +
                "  --humann_protein_db <path>\n" +
                "  (and ensure --metaphlan_db is also provided)"
            )
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Subworkflow 3: ASSEMBLY - Metagenome assembly using MEGAHIT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

        ASSEMBLY(
            FASTP.out.reads
        )
        ch_versions = ch_versions.concat(ASSEMBLY.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nion_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
