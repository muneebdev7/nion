/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nion_pipeline'

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { METAPHLAN3_METAPHLAN3  } from '../modules/nf-core/metaphlan3/metaphlan3/main'
include { METAPHLAN3_MERGEMETAPHLANTABLES } from '../modules/nf-core/metaphlan3/mergemetaphlantables/main'
include { HUMANN                 } from '../modules/local/humann/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NION {

    take:
    ch_raw_short_reads // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_raw_short_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
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


    //
    // MODULE: Run MetaPhlAn3 for taxonomic classification and abundance estimation
    //
    if (params.metaphlan_db) {
        def ch_metaphlan_db
        ch_metaphlan_db = channel.fromPath(params.metaphlan_db, checkIfExists: true).first()

        METAPHLAN3_METAPHLAN3(
            FASTP.out.reads.map { meta, reads -> [meta, reads] },
            ch_metaphlan_db
        )

        ch_versions = ch_versions.mix(METAPHLAN3_METAPHLAN3.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(METAPHLAN3_METAPHLAN3.out.profile.map { _meta, profile -> profile })
    
        METAPHLAN3_MERGEMETAPHLANTABLES (
            METAPHLAN3_METAPHLAN3.out.profile.collect{ _meta, profile -> profile }.map{ profiles -> [[id:'merged'], profiles]}
        )
        ch_versions = ch_versions.mix(METAPHLAN3_MERGEMETAPHLANTABLES.out.versions)
    } else {
        log.warn("MetaPhlAn is disabled: provide --metaphlan_db.")
    }

    //
    // MODULE: Run Humann for functional annotation
    //
    if (params.humann_nucleotide_db && params.humann_protein_db) {
        def ch_humann_nucleotide_db = channel.fromPath(params.humann_nucleotide_db, checkIfExists: true).first()
        def ch_humann_protein_db = channel.fromPath(params.humann_protein_db, checkIfExists: true).first()
        def ch_humann_input = FASTP.out.reads
            .join(METAPHLAN3_METAPHLAN3.out.profile)
            .map { meta, reads, profile -> [meta, reads, profile] }
        
        HUMANN(
            ch_humann_input,
            ch_humann_nucleotide_db,
            ch_humann_protein_db
        )
        ch_versions = ch_versions.mix(HUMANN.out.versions)
    } else {
        log.warn("HUMAnN is disabled: provide --humann_nucleotide_db and --humann_protein_db.")
    }

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
