//
// Non-human contamination screening of BAM/CRAM files with Kraken2 (+ optional Bracken)
//
//   SAMTOOLS_FASTQ   per sample; by default only read pairs with both mates unmapped (-f 12, via ext.args)
//   KRAKEN2_KRAKEN2  per sample, paired-end classification
//   BRACKEN_BRACKEN  per sample, optional abundance re-estimation
//

include { SAMTOOLS_FASTQ  } from '../../../modules/nf-core/samtools/fastq/main'
include { KRAKEN2_KRAKEN2 } from '../../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN } from '../../../modules/nf-core/bracken/bracken/main'

workflow BAM_SCREEN_KRAKEN2 {

    take:
    ch_bam        // channel: [ val(meta), path(bam), path(bai) ]
    ch_kraken2_db // channel: path(db)
    ch_bracken_db // channel: path(db)
    run_bracken   // boolean

    main:

    SAMTOOLS_FASTQ (
        ch_bam.map { meta, bam, _bai -> [ meta + [ single_end: false ], bam ] },
        false
    )

    KRAKEN2_KRAKEN2 ( SAMTOOLS_FASTQ.out.fastq, ch_kraken2_db, false, false )

    def ch_multiqc_files = KRAKEN2_KRAKEN2.out.report.map { _meta, report -> report }
    def ch_bracken = channel.empty()
    if (run_bracken) {
        BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.report, ch_bracken_db )
        ch_bracken = BRACKEN_BRACKEN.out.reports
        ch_multiqc_files = ch_multiqc_files.mix(BRACKEN_BRACKEN.out.reports.map { _meta, report -> report })
    }

    emit:
    kraken2_report = KRAKEN2_KRAKEN2.out.report // channel: [ val(meta), path(report) ]
    bracken_report = ch_bracken                  // channel: [ val(meta), path(tsv) ]
    multiqc_files  = ch_multiqc_files            // channel: path(file)
}
