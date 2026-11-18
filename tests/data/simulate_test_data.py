#!/usr/bin/env python3
"""
Deterministic synthetic test data for nf-core/contaminationqc.

Creates a 1 Mb synthetic reference ('chr21'), a somalier-style sites VCF (INFO/AF) and aligned,
coordinate-sorted, indexed BAMs with a known truth:

  sample  individual  SM tag   truth
  P1_N    P1          P1_N     clean
  P1_T    P1          P1_T     clean (same individual as P1_N -> longitudinal 'ok')
  P2_A    P2          SAMPLE   clean (SM tag shared with P6_A on purpose)
  P3_A    P3          P3_A     5 % contaminated by P2           -> source_identified (P2_A)
  P4_A    P4          P4_A     4 % contaminated by P5 (absent)  -> contaminated_source_not_in_cohort
                               + sarscov2 read pairs (unmapped) -> kraken2 hit
  P6_A    P6          SAMPLE   clean
  P5_X    P5          P5_X     background-only sample (.somalier written by the generator)

Each sample additionally carries random unmapped read pairs (unclassified by kraken2).

Requires numpy and docker (samtools / somalier containers); sarscov2 genome is downloaded from
nf-core/test-datasets. Usage: python simulate_test_data.py [outdir]
"""
import gzip
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np

SAMTOOLS = "community.wave.seqera.io/library/htslib_samtools:1.25--704fb80d75745bd4"
SOMALIER = "community.wave.seqera.io/library/somalier:0.3.6--127e0e29751f4940"
SARSCOV2 = ("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/"
            "genome/genome.fasta")

CONTIG, LENGTH, N_SITES, SPACING, FIRST = "chr21", 1_000_000, 400, 2_000, 1_000
READ_LEN, DEPTH, ERR = 50, 100, 0.001
rng = np.random.default_rng(20261117)
BASES = np.array(list("ACGT"))


def docker(image, cmd, wd):
    subprocess.run(["docker", "run", "--rm", "--platform", "linux/amd64",
                    "-v", f"{wd}:/data", "-w", "/data", image, "bash", "-c", cmd], check=True)


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent).resolve()
    out.mkdir(parents=True, exist_ok=True)

    # ---- reference + sites
    ref = BASES[rng.integers(0, 4, LENGTH)]
    with open(out / "genome.fasta", "w") as fh:
        fh.write(f">{CONTIG}\n")
        s = "".join(ref)
        for i in range(0, LENGTH, 60):
            fh.write(s[i:i + 60] + "\n")
    pos = FIRST + SPACING * np.arange(N_SITES)  # 1-based
    refb = ref[pos - 1]
    altb = np.array([rng.choice([b for b in "ACGT" if b != r]) for r in refb])
    af = np.round(rng.uniform(0.2, 0.8, N_SITES), 4)
    with gzip.open(out / "sites.vcf.gz", "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={CONTIG},length={LENGTH}>\n")
        fh.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for p, r, a, f in zip(pos, refb, altb, af):
            fh.write(f"{CONTIG}\t{p}\t.\t{r}\t{a}\t.\t.\tAF={f}\n")
    # panel BED covering the first 80 % of sites (60 bp windows); tests restrict-sites
    with open(out / "panel.bed", "w") as fh:
        for p in pos[: int(0.8 * N_SITES)]:
            fh.write(f"{CONTIG}\t{p - 31}\t{p + 29}\n")

    # ---- individuals (HWE genotypes, alt allele dosage 0/1/2)
    inds = {k: rng.binomial(2, af) for k in ("P1", "P2", "P3", "P4", "P5", "P6")}
    samples = [  # sample, individual, SM, contaminant individual, alpha, sarscov2 pairs
        ("P1_N", "P1", "P1_N", None, 0.0, 0),
        ("P1_T", "P1", "P1_T", None, 0.0, 0),
        ("P2_A", "P2", "SAMPLE", None, 0.0, 0),
        ("P3_A", "P3", "P3_A", "P2", 0.05, 0),
        ("P4_A", "P4", "P4_A", "P5", 0.04, 400),
        ("P6_A", "P6", "SAMPLE", None, 0.0, 0),
        ("P5_X", "P5", "P5_X", None, 0.0, 0),
    ]
    cov = urllib.request.urlopen(SARSCOV2).read().decode().split("\n", 1)[1].replace("\n", "")
    ref_str = "".join(ref)
    rev = str.maketrans("ACGT", "TGCA")

    for name, ind, sm, cont, alpha, n_cov in samples:
        g = inds[ind] / 2.0
        gc = inds[cont] / 2.0 if cont else np.zeros(N_SITES)
        p_alt = (1 - alpha) * g + alpha * gc
        lines = []
        rid = 0
        for i, p in enumerate(pos):
            n = rng.poisson(DEPTH)
            starts = rng.integers(p - READ_LEN + 1, p + 1, n)  # 1-based leftmost
            alt = rng.random(n) < p_alt[i]
            for st, is_alt in zip(starts, alt):
                seq = list(ref_str[st - 1: st - 1 + READ_LEN])
                if is_alt:
                    seq[p - st] = altb[i]
                errs = np.nonzero(rng.random(READ_LEN) < ERR)[0]
                for e in errs:
                    seq[e] = rng.choice([b for b in "ACGT" if b != seq[e]])
                flag = 16 if rng.random() < 0.5 else 0
                rid += 1
                lines.append((st, f"r{rid}\t{flag}\t{CONTIG}\t{st}\t60\t{READ_LEN}M\t*\t0\t0\t{''.join(seq)}\t"
                                  f"{'I' * READ_LEN}\tRG:Z:rg1"))
        lines.sort(key=lambda x: x[0])
        body = [l for _, l in lines]

        def pair(prefix, k, s1, s2):
            q = "I" * len(s1)
            return [f"{prefix}{k}\t77\t*\t0\t0\t*\t*\t0\t0\t{s1}\t{q}\tRG:Z:rg1",
                    f"{prefix}{k}\t141\t*\t0\t0\t*\t*\t0\t0\t{s2}\t{q}\tRG:Z:rg1"]

        for k in range(100):  # random, unclassifiable unmapped pairs
            s1, s2 = "".join(BASES[rng.integers(0, 4, 100)]), "".join(BASES[rng.integers(0, 4, 100)])
            body += pair("u", k, s1, s2)
        for k in range(n_cov):  # sarscov2 fragments (300 bp inserts, 100 bp mates)
            st = int(rng.integers(0, len(cov) - 300))
            frag = cov[st: st + 300]
            body += pair("v", k, frag[:100], frag[-100:].translate(rev)[::-1])
        header = (f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{CONTIG}\tLN:{LENGTH}\n"
                  f"@RG\tID:rg1\tSM:{sm}\tPL:ILLUMINA\n")
        (out / f"{name}.sam").write_text(header + "\n".join(body) + "\n")
        docker(SAMTOOLS, f"samtools sort -o {name}.bam {name}.sam && samtools index {name}.bam && rm {name}.sam",
               out)
        print(f"wrote {name}.bam  (individual {ind}, SM {sm}, contaminant {cont}, alpha {alpha})")

    # background sample: only its .somalier is part of the test data
    docker(SAMTOOLS, "samtools faidx genome.fasta", out)
    docker(SOMALIER, "somalier extract -d background -s sites.vcf.gz -f genome.fasta "
                     "--sample-prefix P5_X__ P5_X.bam", out)
    for f in ("P5_X.bam", "P5_X.bam.bai"):
        (out / f).unlink()

    with open(out / "samplesheet.csv", "w") as fh:
        fh.write("sample,bam,bai,individual,batch\n")
        for name, ind, *_ in samples[:-1]:
            fh.write(f"{name},tests/data/{name}.bam,tests/data/{name}.bam.bai,{ind},\n")
    with open(out / "background.csv", "w") as fh:
        fh.write("sample,somalier,individual\n")
        fh.write("P5_X,tests/data/background/P5_X.somalier,P5\n")
    with open(out / "baseline_samples.tsv", "w") as fh:  # sample sheet for the original script
        fh.write("bam\tindividual\n")
        for name, ind, *_ in samples[:-1]:
            fh.write(f"{name}.bam\t{ind}\n")


if __name__ == "__main__":
    main()
