# nf-core/contaminationqc: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/contaminationqc/usage](https://nf-co.re/contaminationqc/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/contaminationqc takes aligned human BAM/CRAM files and answers four questions per sample:

1. **Is it contaminated with human DNA, and how much?** (somalier CHARR, mixture-model fraction `alpha_unknown_source`)
2. **Whose DNA is it?** For every receiver R and candidate C, at the sites where R is homozygous-like it compares
   `R = (1-a)·R_self + a·C` (C's genotypes) with `R = (1-a)·R_self + a·random person` (Hardy-Weinberg on the sites' AF).
   `LLR = max_a LL_source - max_a LL_unknown`. LLR >> 0 means C explains R's foreign alleles better than a random person.
3. **Is it who it claims to be?** Hom-hom discordance is robust to a few % contamination, so samples of the same individual
   are recognised even when one of them is contaminated (swaps, longitudinal samples).
4. **Is there non-human DNA?** (optional, Kraken2 on unmapped read pairs)

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with a header row as shown in the example below.

```bash
--input '[path to samplesheet file]'
```

| Column       | Description                                                                                                                                                       |
| ------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`     | Unique sample ID (letters, digits, `.`, `_`, `-`). Used in all outputs.                                                                                           |
| `bam`        | BAM or CRAM file aligned to `--fasta`.                                                                                                                            |
| `bai`        | Index (`.bai`, `.csi`, `.crai`).                                                                                                                                  |
| `individual` | Optional patient/individual ID. Samples of the same individual are expected to be genotype-identical (`longitudinal.pairs.tsv`) and serve as each other's anchor. |
| `batch`      | Optional batch (e.g. sequencing run). Relatedness and source attribution run within a batch; identity checks run across batches.                                  |

```csv title="samplesheet.csv"
sample,bam,bai,individual,batch
PATIENT1_TUMOR,/data/PATIENT1_TUMOR.bam,/data/PATIENT1_TUMOR.bam.bai,PATIENT1,RUN01
PATIENT1_NORMAL,/data/PATIENT1_NORMAL.bam,/data/PATIENT1_NORMAL.bam.bai,PATIENT1,RUN01
PATIENT2_TUMOR,/data/PATIENT2_TUMOR.cram,/data/PATIENT2_TUMOR.cram.crai,PATIENT2,RUN02
```

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

To create a samplesheet from S3 prefixes or local directories, use the helper script shipped with the pipeline (needs the AWS CLI for `s3://`; each location is listed once, indexes are matched from the listing, files without index are skipped):

```bash
bin/bam_dir_to_samplesheet.py s3://bucket/runs/run42/ s3://bucket/runs/run43/ \
    --batch-regex 'runs/([^/]+)/' \
    --individual-regex '^(P[0-9]+)_' \
    --exclude 'Logs_Intermediates' \
    -o samplesheet.csv
```

`--batch-level 1` uses the BAM's parent folder as batch, `--prefix-parent` makes duplicated file names unique, `--profile`/`--region`/`--no-sign-request` are passed to the AWS CLI. See `bin/bam_dir_to_samplesheet.py --help`.

The `SM` tag of the BAM read group is not used as the sample name: somalier stores `<sample>__<SM>` internally so that duplicated `SM` tags (common with vendor pipelines) cannot collide.

### Background samples

`--somalier_background` takes a CSV of previously extracted `.somalier` files (e.g. `somalier/extract/` of earlier runs) that act as additional candidate sources and identity references but are never reported as samples themselves:

```csv title="background.csv"
sample,somalier,individual
RUN0815_PATIENT9,/archive/run0815/somalier/extract/RUN0815_PATIENT9.somalier,PATIENT9
```

They must have been extracted with the same sites (after panel restriction) as the current run, and their IDs and file names must not clash with the current samples.

## Reference data

- `--fasta` (or `--genome` with iGenomes): the reference the BAMs are aligned to; `--fasta_fai` is generated when not given.
- `--sites`: somalier sites VCF with `INFO/AF` in the same build and contig naming. The pipeline fails early if no contig name is shared with the FASTA index.
- `--panel_bed` / `--bed_pad`: restrict sites to (padded) panel targets. Strongly recommended for targeted panels; with fewer than ~300 remaining sites consider building a panel-specific sites file with `somalier find-sites`.

## Interpreting calls

| `call`                              | Meaning                                                                                                     |
| ----------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `clean`                             | No evidence of contamination above `--contam_threshold`.                                                    |
| `source_identified`                 | A cohort/background sample explains the foreign alleles (`LLR >= --min_llr`, margin `>= --min_margin`).     |
| `source_ambiguous`                  | Several individuals explain the contamination similarly well (e.g. related candidates).                     |
| `contaminated_source_not_in_cohort` | Contamination `>= --contam_threshold` but no candidate explains it: the source was not sequenced with it.   |
| `…;low_sites`                       | Fewer than 50 homozygous-like sites were usable; interpret with care (low coverage or too few panel sites). |

Samples of the same genotype group (same individual) are never reported as each other's contamination source.

## Large cohorts

The pipeline was designed for thousands of samples per run:

- `somalier extract` runs per sample; relatedness and source attribution run per `batch`; attribution is additionally split into chunks of `--attribution_chunk_size` receivers (default 250).
- The anchored `somalier contamination` MLE is computed only within genotype groups by default. `--somalier_anchored_all_pairs` restores the all-vs-all behaviour of the original script (quadratic in batch size).
- `--report_top_candidates 25` keeps `attribution.pairs.tsv.gz` small (it otherwise has N×(N-1) rows per batch).
- Benchmark (synthetic, 3,000 samples in one batch, 1,000 sites, 4 CPUs): identity 2 s, attribution 20 s per 250-receiver chunk, `somalier relate` 10 s, report 40 s with 7.7 GB RAM. All 30 simulated contaminations (1–5 %) were attributed to the correct individual, with no false positives among 2,970 clean samples. See `tests/benchmark/simulate_cohort.py`.
- On small panels the somalier relatedness estimate is noisy; raise `--relatedness_flag_threshold` (e.g. to 0.35) to avoid thousands of `elevated_relatedness` flags.

## Kraken2 screening

With `--run_kraken2 --kraken2_db <dir or .tar.gz>`, read pairs with both mates unmapped (`--kraken2_reads unmapped`, default) or all primary reads (`--kraken2_reads all`) are extracted with `samtools fastq` and classified with Kraken2; `--run_bracken` adds Bracken abundance estimates. Paired-end BAMs are expected. For CRAM input with `--kraken2_reads all`, the reference must be resolvable by samtools (e.g. `REF_PATH`/`REF_CACHE`), as `samtools fastq` is not given the FASTA.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/contaminationqc --input ./samplesheet.csv --outdir ./results --fasta genome.fa --sites sites.vcf.gz -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/contaminationqc -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/contaminationqc
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/contaminationqc releases page](https://github.com/nf-core/contaminationqc/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
