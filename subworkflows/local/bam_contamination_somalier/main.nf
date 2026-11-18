//
// Contamination, contamination-source attribution and sample identity with somalier
//
//   SOMALIER_EXTRACT            per sample
//   SOMALIER_RELATE             per batch
//   SOMALIERCONTAM_IDENTITY     global (cohort + background) -> genotype groups, identity pairs
//   SOMALIER_CONTAMINATION      per genotype group (default) or per batch (all ordered pairs)
//   SOMALIERCONTAM_ATTRIBUTE    per batch x receiver shard (candidates: batch cohort + background)
//   SOMALIERCONTAM_REPORT       global
//

include { SOMALIER_EXTRACT             } from '../../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE              } from '../../../modules/nf-core/somalier/relate/main'
include { SOMALIER_CONTAMINATION       } from '../../../modules/local/somalier/contamination/main'
include { SOMALIERCONTAM_IDENTITY      } from '../../../modules/local/somaliercontam/identity/main'
include { SOMALIERCONTAM_ATTRIBUTE     } from '../../../modules/local/somaliercontam/attribute/main'
include { SOMALIERCONTAM_REPORT        } from '../../../modules/local/somaliercontam/report/main'

// Header of the sample table read by somalier_contam.py
def sampleTableHeader() {
    return "sample\tfile\tindividual\tbatch\trole\tinput"
}

// Sample name stored inside a .somalier file (format v2: version byte, name length byte, name)
def somalierSampleName(somalier) {
    return somalier.withInputStream { stream ->
        stream.read()
        def n = stream.read()
        return new String(stream.readNBytes(n), 'UTF-8')
    }
}

// One row of the sample table read by somalier_contam.py. Files are staged as somalier/<name>.
def sampleTableRow(meta, somalier, input_name) {
    return [ meta.id, "somalier/${somalier.name}", meta.individual ?: '', meta.batch ?: '', meta.role, input_name ?: '' ].join('\t')
}

workflow BAM_CONTAMINATION_SOMALIER {

    take:
    ch_bam                 // channel: [ val(meta), path(bam), path(bai) ]  meta: id, individual, batch
    ch_background          // channel: [ val(meta), path(somalier) ]        meta: id, individual, role
    ch_fasta               // channel: [ val(meta), path(fasta) ]
    ch_fai                 // channel: [ val(meta), path(fai) ]
    ch_sites               // channel: [ val(meta), path(vcf) ]
    attribution_chunk_size // integer: receivers per attribution task
    anchored_all_pairs     // boolean: run somalier contamination on all ordered pairs per batch

    main:

    //
    // MODULE: somalier extract (sample names are prefixed with the pipeline sample ID via ext.args)
    //
    SOMALIER_EXTRACT ( ch_bam, ch_fasta, ch_fai, ch_sites )

    def ch_cohort = SOMALIER_EXTRACT.out.extract
        .join(ch_bam.map { meta, bam, _bai -> [ meta, bam.name ] })
        .map { meta, somalier, input_name -> [ meta + [ role: 'cohort' ], somalier, input_name ] }
    def ch_bg = ch_background
        .map { meta, somalier -> [ meta + [ batch: '', role: 'background' ], somalier, '' ] }
    def ch_all = ch_cohort.mix(ch_bg)

    //
    // Global sample table + files (cohort and background)
    //
    def ch_all_table = ch_all
        .map { meta, somalier, input_name -> sampleTableRow(meta, somalier, input_name) }
        .collectFile(name: 'samples.tsv', seed: sampleTableHeader() + '\n', newLine: true, sort: true)
    def ch_all_files = ch_all.map { _meta, somalier, _input -> somalier }.collect(sort: { a, b -> a.name <=> b.name })

    //
    // MODULE: identity (hom-hom discordance, genotype groups) across all samples and batches
    //
    SOMALIERCONTAM_IDENTITY (
        ch_all_files.map { files -> [ files ] }.combine(ch_all_table).map { files, table ->
            [ [ id: 'identity' ], files, table ]
        },
        ch_sites
    )

    //
    // Batches of cohort samples: [ batch, [metas], [somalier], [input names] ]
    //
    def ch_batches = ch_cohort
        .map { meta, somalier, input_name -> [ meta.batch ?: 'cohort', meta, somalier, input_name ] }
        .groupTuple()
        // deterministic sample order within a batch (stable outputs and -resume caching)
        .map { batch, metas, somaliers, inputs ->
            def rows = [ metas, somaliers, inputs ].transpose().sort { row -> row[0].id }
            [ batch, rows.collect { row -> row[0] }, rows.collect { row -> row[1] }, rows.collect { row -> row[2] } ]
        }

    //
    // MODULE: somalier relate per batch, with expected groups (individuals with > 1 sample)
    //
    def ch_groups_files = ch_batches
        .flatMap { batch, metas, somaliers, _inputs ->
            [ metas, somaliers ].transpose()
                .findAll { meta, _somalier -> meta.individual }
                .groupBy { meta, _somalier -> meta.individual }
                .values()
                .findAll { members -> members.size() > 1 }
                .collect { members -> [ "${batch}.groups.txt", members.collect { _meta, somalier -> somalierSampleName(somalier) }.join(',') + '\n' ] }
        }
        .collectFile()
        .map { groups -> [ groups.name - '.groups.txt', groups ] }

    def ch_relate_in = ch_batches
        .join(ch_groups_files, remainder: true)
        .filter { row -> row[1] != null }
        .combine(ch_sites)
        .multiMap { batch, _metas, somaliers, _inputs, groups, sites_meta, sites ->
            extract: [ [ id: batch ], somaliers, [] ]
            groups:  groups ?: []
            sites:   [ sites_meta, sites ]
        }
    SOMALIER_RELATE ( ch_relate_in.extract, ch_relate_in.groups, ch_relate_in.sites )

    //
    // MODULE: somalier contamination (anchored MLEs)
    //   default: within genotype groups (same individual, incl. background), all pairs otherwise
    //
    def ch_contamination_in
    if (anchored_all_pairs) {
        ch_contamination_in = ch_batches.map { batch, _metas, somaliers, _inputs -> [ [ id: batch ], somaliers ] }
    } else {
        def ch_somalier_by_id = ch_all.map { meta, somalier, _input -> [ meta.id, somalier ] }
        ch_contamination_in = SOMALIERCONTAM_IDENTITY.out.groups
            .flatMap { _meta, groups -> groups.splitCsv(header: true, sep: '\t').collect { row -> [ row.sample, row.genotype_group ] } }
            .join(ch_somalier_by_id)
            .map { _id, group, somalier -> [ group, somalier ] }
            .groupTuple()
            .filter { _group, somaliers -> somaliers.size() > 1 }
            .map { group, somaliers -> [ [ id: group ], somaliers.sort { somalier -> somalier.name } ] }
    }
    SOMALIER_CONTAMINATION ( ch_contamination_in, ch_sites )

    //
    // MODULE: source attribution per batch x receiver shard
    //
    def ch_batch_rows = ch_batches
        .combine(ch_bg.toSortedList { a, b -> a[0].id <=> b[0].id }.map { bg -> [ bg ] })
        .map { batch, metas, somaliers, inputs, bg -> [ batch, metas.size(), [ metas, somaliers, inputs ].transpose() + bg ] }

    def ch_batch_tables = ch_batch_rows
        .map { batch, _n, rows ->
            [ "${batch}.attribution_samples.tsv", ([ sampleTableHeader() ] + rows.collect { meta, somalier, input_name -> sampleTableRow(meta, somalier, input_name) }).join('\n') + '\n' ]
        }
        .collectFile()
        .map { table -> [ table.name - '.attribution_samples.tsv', table ] }

    def ch_attribution_in = ch_batch_rows
        .join(ch_batch_tables)
        .flatMap { batch, n_cohort, rows, table ->
            def n_shards = Math.max(1, Math.ceil(n_cohort / attribution_chunk_size) as int)
            def files = rows.collect { _meta, somalier, _input -> somalier }
            (0..<n_shards).collect { k ->
                [ [ id: "${batch}.shard${k}", batch: batch, shard_index: k, n_shards: n_shards ], files, table ]
            }
        }
    SOMALIERCONTAM_ATTRIBUTE ( ch_attribution_in, ch_sites )

    //
    // MODULE: report (global)
    //
    SOMALIERCONTAM_REPORT (
        ch_all_table
            .combine(SOMALIERCONTAM_IDENTITY.out.groups.map { _meta, groups -> groups })
            .combine(SOMALIERCONTAM_IDENTITY.out.pairs.map { _meta, pairs -> pairs })
            .map { table, groups, pairs -> [ [ id: 'contaminationqc' ], table, groups, pairs ] },
        SOMALIERCONTAM_ATTRIBUTE.out.pairs.map { _meta, pairs -> pairs }.collect(sort: true).map { pairs -> [ pairs ] }
            .combine(SOMALIERCONTAM_ATTRIBUTE.out.self.map { _meta, self -> self }.collect(sort: true).map { self -> [ self ] })
            .map { pairs, self -> [ [ id: 'attribution' ], pairs, self ] },
        SOMALIER_RELATE.out.samples_tsv.map { _meta, tsv -> tsv }.collect(sort: true).map { tsv -> [ tsv ] }
            .combine(SOMALIER_RELATE.out.pairs_tsv.map { _meta, tsv -> tsv }.collect(sort: true).map { tsv -> [ tsv ] })
            .map { samples, pairs -> [ [ id: 'relate' ], samples, pairs ] },
        SOMALIER_CONTAMINATION.out.pairs_tsv.map { _meta, tsv -> tsv }.collect(sort: true).ifEmpty([])
            .map { pairs -> [ [ id: 'contamination' ], pairs ] }
    )

    def ch_multiqc_files = SOMALIERCONTAM_REPORT.out.multiqc.map { _meta, files -> files }.flatten()
        .mix(SOMALIER_RELATE.out.samples_tsv.map { _meta, tsv -> tsv })
        .mix(SOMALIER_RELATE.out.pairs_tsv.map { _meta, tsv -> tsv })

    emit:
    extract           = SOMALIER_EXTRACT.out.extract               // channel: [ val(meta), path(somalier) ]
    summary           = SOMALIERCONTAM_REPORT.out.summary          // channel: [ val(meta), path(tsv) ]
    attribution_pairs = SOMALIERCONTAM_REPORT.out.attribution_pairs // channel: [ val(meta), path(tsv.gz) ]
    longitudinal      = SOMALIERCONTAM_REPORT.out.longitudinal     // channel: [ val(meta), path(tsv) ]
    flags             = SOMALIERCONTAM_REPORT.out.flags            // channel: [ val(meta), path(tsv) ]
    relate_html       = SOMALIER_RELATE.out.html                   // channel: [ val(meta), path(html) ]
    multiqc_files     = ch_multiqc_files                           // channel: path(file)
}
