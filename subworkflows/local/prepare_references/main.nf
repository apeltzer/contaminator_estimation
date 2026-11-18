//
// Prepare global reference inputs: FASTA index, (panel-restricted) somalier sites, Kraken2 database
//

include { SAMTOOLS_FAIDX               } from '../../../modules/nf-core/samtools/faidx/main'
include { UNTAR as UNTAR_KRAKEN2_DB    } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_BRACKEN_DB    } from '../../../modules/nf-core/untar/main'
include { SOMALIERCONTAM_RESTRICTSITES } from '../../../modules/local/somaliercontam/restrictsites/main'
include { checkSitesContigs            } from '../utils_nfcore_contaminationqc_pipeline'

workflow PREPARE_REFERENCES {

    take:
    fasta       //  string: path to reference FASTA
    fasta_fai   //  string: path to FASTA index (optional)
    sites       //  string: path to somalier sites VCF with INFO/AF
    panel_bed   //  string: path to panel target BED (optional)
    kraken2_db  //  string: path to Kraken2 database directory or .tar.gz (optional)
    bracken_db  //  string: path to Bracken-enabled database directory or .tar.gz (optional)

    main:

    def ch_fasta = channel.value([ [ id: file(fasta).baseName ], file(fasta, checkIfExists: true) ])

    //
    // FASTA index
    //
    def ch_fai
    if (fasta_fai) {
        ch_fai = channel.value([ [ id: file(fasta).baseName ], file(fasta_fai, checkIfExists: true) ])
    } else {
        SAMTOOLS_FAIDX ( ch_fasta.map { meta, fa -> [ meta, fa, [] ] }, false )
        ch_fai = SAMTOOLS_FAIDX.out.fai.first()
    }

    //
    // Sites: sanity-check contig naming against the FASTA, optionally restrict to panel targets
    //
    def sites_file = file(sites, checkIfExists: true)
    def ch_sites_raw = ch_fai
        .map { _meta, fai -> checkSitesContigs(sites_file, fai) }
        .map { _ok -> [ [ id: sites_file.name.replaceAll(/\.vcf(\.b?gz)?$/, '') ], sites_file ] }

    def ch_sites
    if (panel_bed) {
        SOMALIERCONTAM_RESTRICTSITES (
            ch_sites_raw,
            channel.value([ [ id: file(panel_bed).baseName ], file(panel_bed, checkIfExists: true) ])
        )
        ch_sites = SOMALIERCONTAM_RESTRICTSITES.out.vcf.first()
    } else {
        ch_sites = ch_sites_raw.first()
    }

    //
    // Kraken2 / Bracken databases (directory or tar.gz archive)
    //
    def ch_kraken2_db = channel.empty()
    if (kraken2_db) {
        if (kraken2_db.toString().endsWith('.tar.gz') || kraken2_db.toString().endsWith('.tgz')) {
            UNTAR_KRAKEN2_DB ( channel.of([ [ id: 'kraken2_db' ], file(kraken2_db, checkIfExists: true) ]) )
            ch_kraken2_db = UNTAR_KRAKEN2_DB.out.untar.map { _meta, db -> db }.first()
        } else {
            ch_kraken2_db = channel.value(file(kraken2_db, checkIfExists: true))
        }
    }
    def ch_bracken_db = ch_kraken2_db
    if (bracken_db && bracken_db != kraken2_db) {
        if (bracken_db.toString().endsWith('.tar.gz') || bracken_db.toString().endsWith('.tgz')) {
            UNTAR_BRACKEN_DB ( channel.of([ [ id: 'bracken_db' ], file(bracken_db, checkIfExists: true) ]) )
            ch_bracken_db = UNTAR_BRACKEN_DB.out.untar.map { _meta, db -> db }.first()
        } else {
            ch_bracken_db = channel.value(file(bracken_db, checkIfExists: true))
        }
    }

    emit:
    fasta      = ch_fasta      // channel: [ val(meta), path(fasta) ]
    fai        = ch_fai        // channel: [ val(meta), path(fai) ]
    sites      = ch_sites      // channel: [ val(meta), path(vcf) ]
    kraken2_db = ch_kraken2_db // channel: path(db)
    bracken_db = ch_bracken_db // channel: path(db)
}
