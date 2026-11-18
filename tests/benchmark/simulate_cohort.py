#!/usr/bin/env python3
"""
Simulate a large cohort of .somalier files (format v2) for scale benchmarks of the
somalier_contam.py steps / pipeline, without BAMs.

  simulate_cohort.py <outdir> [n_samples=3000] [n_sites=1000] [depth=400] [n_contaminated=30]

Each sample is a distinct individual (HWE genotypes on the sites' AF), except that 10 % of individuals
(every 9th individual) has a second sample of the same individual (longitudinal pairs). n_contaminated samples receive 1-5 %
DNA from another cohort sample (truth written to truth.tsv). Writes sites.vcf, samples.tsv (sample
table for somalier_contam.py) and samplesheet-like background.csv.
"""
import struct
import sys
from pathlib import Path

import numpy as np


def write_somalier(path, name, counts):
    with open(path, "wb") as fh:
        fh.write(struct.pack("<BB", 2, len(name)) + name.encode())
        fh.write(struct.pack("<HHH", len(counts), 0, 0))
        fh.write(counts.astype("<u4").tobytes())


def main():
    out = Path(sys.argv[1])
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    s = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    depth = int(sys.argv[4]) if len(sys.argv) > 4 else 400
    n_cont = int(sys.argv[5]) if len(sys.argv) > 5 else 30
    rng = np.random.default_rng(1)
    (out / "somalier").mkdir(parents=True, exist_ok=True)
    af = rng.uniform(0.15, 0.85, s)
    # sites: alleles A/C so REF < ALT (no AF flip), sorted single contig
    with open(out / "sites.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=250000000>\n"
                 '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for i, f in enumerate(af):
            fh.write(f"chr1\t{1000 + i * 5000}\t.\tA\tC\t.\t.\tAF={f:.4f}\n")
    n_ind = n - n // 10
    geno = rng.binomial(2, af, size=(n_ind, s)) / 2.0
    ind_of = list(range(n_ind)) + [(i * 9) % n_ind for i in range(n // 10)]
    cont_src = {}
    receivers = rng.choice(n, n_cont, replace=False)
    for r in receivers:
        src = int(rng.integers(0, n))
        while ind_of[src] == ind_of[r]:
            src = int(rng.integers(0, n))
        cont_src[int(r)] = (src, float(rng.uniform(0.01, 0.05)))
    with open(out / "samples.tsv", "w") as tab, open(out / "truth.tsv", "w") as tr:
        tab.write("sample\tfile\tindividual\tbatch\trole\tinput\n")
        tr.write("sample\tsource\talpha\n")
        for k in range(n):
            name = f"S{k:05d}"
            p = geno[ind_of[k]].copy()
            if k in cont_src:
                src, a = cont_src[k]
                p = (1 - a) * p + a * geno[ind_of[src]]
                tr.write(f"{name}\tS{src:05d}\t{a:.4f}\n")
            dp = rng.poisson(depth, s)
            alt = rng.binomial(dp, np.clip(p * 0.998 + 0.001, 0, 1))
            write_somalier(out / "somalier" / f"{name}.somalier", f"{name}__{name}",
                           np.stack([dp - alt, alt, np.zeros(s, dtype=int)], axis=1))
            tab.write(f"{name}\tsomalier/{name}.somalier\tI{ind_of[k]:05d}\t\tcohort\t{name}.bam\n")
    print(f"wrote {n} samples x {s} sites to {out}")


if __name__ == "__main__":
    main()
