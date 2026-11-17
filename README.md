<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-contaminationqc_logo_dark.png">
    <img alt="nf-core/contaminationqc" src="docs/images/nf-core-contaminationqc_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/contaminationqc)
[![GitHub Actions CI Status](https://github.com/nf-core/contaminationqc/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/contaminationqc/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/contaminationqc/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/contaminationqc/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/contaminationqc/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.2.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.2.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/contaminationqc)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23contaminationqc-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/contaminationqc)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/contaminationqc** is a bioinformatics pipeline that checks aligned human sequencing data (BAM/CRAM, e.g. targeted panels such as TSO500) for **sample contamination**, identifies **which sample in the cohort the contaminating DNA came from**, verifies **sample identity** (swaps, longitudinal same-patient samples) and optionally screens for **non-human contamination** with Kraken2. It scales to thousands of samples per run.

1. Prepare references: FASTA index ([`SAMtools`](https://www.htslib.org/)), restrict somalier sites to panel targets (optional)
2. Extract allele counts at informative sites per sample ([`somalier extract`](https://github.com/brentp/somalier))
3. Relatedness, CNV-robust concordance and CHARR contamination per batch ([`somalier relate`](https://github.com/brentp/somalier))
4. Sample identity: hom-hom discordance and genotype (= individual) groups across all samples and batches (`somalier_contam.py identity`)
5. Anchored contamination MLEs between samples of the same individual ([`somalier contamination`](https://github.com/brentp/somalier))
6. Contamination **source attribution**: mixture model "candidate sample is the contaminant" vs "an unknown person is the contaminant", sharded over receivers (`somalier_contam.py attribute`)
7. Per-sample calls (`clean`, `source_identified`, `source_ambiguous`, `contaminated_source_not_in_cohort`), longitudinal and relatedness flags (`somalier_contam.py report`)
8. Optional: non-human screening of unmapped read pairs ([`SAMtools fastq`](https://www.htslib.org/), [`Kraken2`](https://ccb.jhu.edu/software/kraken2/), [`Bracken`](https://ccb.jhu.edu/software/bracken/))
9. Present QC ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,bam,bai,individual,batch
PATIENT1_TUMOR,/data/PATIENT1_TUMOR.bam,/data/PATIENT1_TUMOR.bam.bai,PATIENT1,RUN01
PATIENT1_NORMAL,/data/PATIENT1_NORMAL.bam,/data/PATIENT1_NORMAL.bam.bai,PATIENT1,RUN01
PATIENT2_TUMOR,/data/PATIENT2_TUMOR.cram,/data/PATIENT2_TUMOR.cram.crai,PATIENT2,RUN02
```

Each row is one BAM/CRAM. `individual` (optional) groups samples of the same patient: they are checked for identity and used as clean anchors. `batch` (optional, e.g. the sequencing run or library-prep plate) limits the search for contamination sources to samples processed together, which keeps the all-vs-all comparison tractable for thousands of samples. Identity checks always run across all batches.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/contaminationqc \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --fasta /refs/hg19.fa \
   --sites /refs/somalier/sites.hg19.vcf.gz \
   --panel_bed /refs/tso500_v2_targets.bed \
   --outdir <OUTDIR>
```

### Screening runs with predefined sites and a Kraken2 database

For routine screening (e.g. a nightly run over all samples of the day), keep the reference resources fixed in a parameter file and only change the samplesheet and output directory:

`screening.yaml`:

```yaml
fasta: /refs/hg19/genome.fa
fasta_fai: /refs/hg19/genome.fa.fai
sites: /refs/somalier/sites.hg19.vcf.gz # predefined somalier sites with INFO/AF, same build + contig names as the BAMs
panel_bed: /refs/tso500_v2_targets.bed # optional: restrict sites to the panel
run_kraken2: true
kraken2_db: /refs/kraken2/k2_pluspf_16gb # directory or .tar.gz of a prebuilt Kraken2 database
kraken2_memory_mapping: true # database on shared storage is memory-mapped instead of loaded per task
kraken2_reads: unmapped # classify read pairs with both mates unmapped
run_bracken: true # optional, needs a Bracken-enabled database (all prebuilt k2 databases are)
somalier_background: /screening/background.csv # optional: .somalier files of earlier runs as extra candidate sources
report_top_candidates: 25 # keep the 25 best candidate sources per sample in attribution.pairs.tsv.gz
relatedness_flag_threshold: 0.35 # fewer noise-level relatedness flags on small panels with many samples
```

```bash
nextflow run nf-core/contaminationqc \
   -profile <docker/singularity/.../institute> \
   -params-file screening.yaml \
   --input samplesheet_2026-11-17.csv \
   --outdir /screening/2026-11-17
```

- **Sites:** use a somalier sites VCF for your genome build that has `INFO/AF` (the contamination models need population allele frequencies), e.g. the sites files from the [somalier releases](https://github.com/brentp/somalier/releases) or a panel-specific set made with `somalier find-sites`. Check with `zgrep -v '^#' sites.vcf.gz | head -1 | grep AF=`. Contig names must match the BAMs and FASTA (`chr1` vs `1`); the pipeline checks this at start-up.
- **Kraken2 database:** any Kraken2 database works; prebuilt ones (Standard, PlusPF, capped 8/16 GB versions) are available from the [Kraken2 index zone](https://benlangmead.github.io/aws-indexes/k2). Extract the archive once and point `--kraken2_db` at the directory so it is not unpacked for every run. With thousands of samples, a capped database plus `--kraken2_memory_mapping` keeps memory per task low.
- **Background:** the `somalier/extract/*.somalier` files of earlier runs can be listed in a background CSV (`sample,somalier,individual`) so that a contaminant processed on a different day can still be identified. They must have been extracted with the same (panel-restricted) sites.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/contaminationqc/usage) and the [parameter documentation](https://nf-co.re/contaminationqc/parameters).

## Pipeline output

The main result is `contamination/contaminationqc.summary.samples.tsv` with one row per sample: contamination estimates (CHARR, mixture model), the most likely source sample with its contamination fraction and log-likelihood ratio, and a call. To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/contaminationqc/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/contaminationqc/output).

## Credits

nf-core/contaminationqc was originally written by Alexander Peltzer.

We thank the following people for their extensive assistance in the development of this pipeline:


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#contaminationqc` channel](https://nfcore.slack.com/channels/contaminationqc) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/contaminationqc for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
