/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                    } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_COVERAGE          } from '../modules/nf-core/samtools/coverage/main'
include { PREPARE_REFERENCES         } from '../subworkflows/local/prepare_references'
include { BAM_CONTAMINATION_SOMALIER } from '../subworkflows/local/bam_contamination_somalier'
include { BAM_SCREEN_KRAKEN2         } from '../subworkflows/local/bam_screen_kraken2'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_contaminationqc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CONTAMINATIONQC {

    take:
    ch_samplesheet // channel: [ val(meta), path(bam), path(bai) ] read in from --input (meta: id, individual, batch, sex)
    ch_background  // channel: [ val(meta), path(somalier) ] read in from --somalier_background
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: reference preparation (FASTA index, panel-restricted sites, Kraken2 DB)
    //
    PREPARE_REFERENCES (
        params.fasta,
        params.fasta_fai,
        params.sites,
        params.panel_bed,
        params.run_kraken2 ? params.kraken2_db : null,
        params.run_bracken ? (params.bracken_db ?: params.kraken2_db) : null
    )

    //
    // MODULE: per-contig read counts (MAPQ-filtered) for the sex check
    //
    def ch_coverage = channel.empty()
    if (!params.skip_sex_check) {
        SAMTOOLS_COVERAGE (
            ch_samplesheet,
            PREPARE_REFERENCES.out.fasta
                .combine(PREPARE_REFERENCES.out.fai)
                .map { meta, fasta, _meta2, fai -> [ meta, fasta, fai ] }
                .first()
        )
        ch_coverage = SAMTOOLS_COVERAGE.out.coverage
    }

    //
    // SUBWORKFLOW: somalier relatedness, contamination, source attribution, identity and sex check
    //
    BAM_CONTAMINATION_SOMALIER (
        ch_samplesheet,
        ch_background,
        PREPARE_REFERENCES.out.fasta,
        PREPARE_REFERENCES.out.fai,
        PREPARE_REFERENCES.out.sites,
        ch_coverage,
        params.attribution_chunk_size,
        params.somalier_anchored_all_pairs
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_CONTAMINATION_SOMALIER.out.multiqc_files)

    //
    // SUBWORKFLOW: optional non-human screening with Kraken2 (+ Bracken)
    //
    if (params.run_kraken2) {
        BAM_SCREEN_KRAKEN2 (
            ch_samplesheet,
            PREPARE_REFERENCES.out.kraken2_db,
            PREPARE_REFERENCES.out.bracken_db,
            params.run_bracken
        )
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SCREEN_KRAKEN2.out.multiqc_files)
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

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'contaminationqc_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'contaminationqc'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    summary        = BAM_CONTAMINATION_SOMALIER.out.summary                       // channel: [ val(meta), path(tsv) ]
    versions       = ch_versions                                                  // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
