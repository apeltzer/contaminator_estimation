# nf-core/contaminationqc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/contaminationqc, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Port of the standalone `somalier_contam.py` (relatedness, contamination, contamination-source attribution, identity) to a modular Nextflow pipeline; the original is kept in `tests/baseline/` and results are verified to be equivalent on the test data.
- Per-sample `somalier extract`, per-batch `somalier relate`, global identity, per-group/per-batch `somalier contamination`, sharded source attribution and a global report.
- `batch` samplesheet column and `--somalier_background` for scaling to thousands of samples and identifying contaminants from earlier runs.
- Optional Kraken2 (+ Bracken) screening of unmapped read pairs.
- Sex check: optional samplesheet column `sex`, inferred sex from MAPQ-filtered chrY reads (nf-core `samtools/coverage`) and X heterozygosity (somalier), `sex_check` column in the summary and MultiQC table; `--skip_sex_check`, `--sex_check_min_mapq`, `--sex_female_max_y_rpm`, `--sex_male_min_y_rpm`. Synthetic test data now include chrX/chrY.
- Metro map (`docs/images/contaminationqc_metro_map.mmd`, rendered with [nf-metro](https://github.com/seqeralabs/nf-metro)) and a detailed Mermaid workflow diagram in the README.
- Local modules `somalier/contamination`, `somaliercontam/{restrictsites,identity,attribute,report}`; patched nf-core `somalier/relate` (sites input) and `somalier/extract` (output named by prefix).

### `Fixed`

### `Dependencies`

### `Deprecated`
