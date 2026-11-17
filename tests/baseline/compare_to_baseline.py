#!/usr/bin/env python3
"""
Compare nf-core/contaminationqc outputs with the original all-in-one somalier_contam.py.

  compare_to_baseline.py <baseline_outdir> <pipeline_outdir> [--all-pairs]

Sample names: the original names samples by SM tag (duplicates re-extracted as '<bam stem>-SM');
the pipeline uses the samplesheet ID. The mapping is taken from the original summary's 'bam' column.
Floats are compared with rtol=1e-6/atol=1e-9 (pipeline writes %.4g, as the original), integers and
calls exactly. The best source is compared at the individual level (best_source_individual) and only where its
LLR is not numerical noise (|LLR| > 1e-6): exact LLR ties between samples of the same individual
(and noise-level ties for clean samples) are broken by an unstable sort in the original and by
sample order in the pipeline.
--all-pairs: pipeline was run with --somalier_anchored_all_pairs, compare somalier_anchored_mle per pair.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FLOAT_TOL = dict(rtol=1e-6, atol=1e-9)


def fail(msg):
    print(f"FAIL: {msg}")
    return 1


def cmp_frame(name, a, b, key, cols):
    errors = 0
    m = a.merge(b, on=key, how="outer", suffixes=("_orig", "_nf"), indicator=True)
    if (m._merge != "both").any():
        errors += fail(f"{name}: rows only in one side:\n{m[m._merge != 'both'][key]}")
    m = m[m._merge == "both"]
    for c in cols:
        x, y = m[f"{c}_orig"], m[f"{c}_nf"]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            ok = np.isclose(x.astype(float), y.astype(float), equal_nan=True, **FLOAT_TOL)
        else:
            ok = x.fillna("").astype(str).values == y.fillna("").astype(str).values
        if not ok.all():
            errors += fail(f"{name}.{c}: {int((~ok).sum())} mismatches, e.g.\n"
                           f"{m.loc[~ok, key + [f'{c}_orig', f'{c}_nf']].head()}")
    print(f"{name}: {len(m)} rows, {len(cols)} columns compared, {'OK' if not errors else 'MISMATCH'}")
    return errors


def main():
    orig, nf = Path(sys.argv[1]), Path(sys.argv[2])
    all_pairs = "--all-pairs" in sys.argv
    o_sum = pd.read_csv(orig / "summary.samples.tsv", sep="\t")
    n_sum = pd.read_csv(next(nf.glob("contamination/*.summary.samples.tsv")), sep="\t")
    name_map = {s: Path(b).stem for s, b in zip(o_sum["sample"], o_sum["bam"])}
    rn = lambda s: s.map(lambda x: name_map.get(x, x))  # noqa: E731
    o_sum["sample"] = rn(o_sum["sample"])
    o_sum["best_source"] = rn(o_sum["best_source"])
    errors = 0

    cols = ["individual", "call", "contamination_charr", "alpha_unknown_source", "p_middling_ab",
            "somalier_anchored_mle_same_individual", "n_hom_sites", "llr_contam_vs_clean", "depth_mean",
            "gt_depth_mean", "n_het", "n_unknown"]
    errors += cmp_frame("summary", o_sum, n_sum, ["sample"], [c for c in cols if c in o_sum and c in n_sum])
    sig = o_sum.set_index("sample").source_llr_vs_unknown.abs() > 1e-6
    # best-source columns only where the best source is not an exact/noise-level tie
    src_cols = ["source_alpha", "source_llr_vs_unknown", "llr_margin", "source_n_informative", "second_best_llr"]
    keep = sig[sig].index
    errors += cmp_frame("summary(best source)", o_sum[o_sum["sample"].isin(keep)], n_sum[n_sum["sample"].isin(keep)],
                        ["sample"], src_cols)
    bs = pd.DataFrame({"orig": o_sum.set_index("sample").best_source_individual[sig],
                       "nf": n_sum.set_index("sample").best_source_individual})
    bs = bs[bs.orig.notna()]
    if (bs.orig != bs.nf).any():
        errors += fail(f"best_source_individual mismatch:\n{bs}")
    else:
        print(f"best_source_individual: {len(bs)} non-trivial calls identical")

    o_p = pd.read_csv(orig / "attribution.pairs.tsv", sep="\t")
    n_p = pd.read_csv(next(nf.glob("contamination/*.attribution.pairs.tsv.gz")), sep="\t")
    for c in ("receiver", "candidate"):
        o_p[c] = rn(o_p[c])
    pcols = ["n_sites", "n_informative", "alpha_source", "alpha_unknown_same_sites", "llr_vs_unknown",
             "llr_vs_clean", "hom_discordance", "n_both_hom", "same_genotype_group", "somalier_relatedness",
             "somalier_concordance", "somalier_ibs0", "somalier_n"] + (["somalier_anchored_mle"] if all_pairs else [])
    errors += cmp_frame("attribution.pairs", o_p, n_p, ["receiver", "candidate"], pcols)

    def unordered(df):
        df = df.copy()
        for c in ("sample_a", "sample_b"):
            df[c] = rn(df[c])
        a, b = df[["sample_a", "sample_b"]].min(axis=1), df[["sample_a", "sample_b"]].max(axis=1)
        df["sample_a"], df["sample_b"] = a, b
        return df

    o_l = unordered(pd.read_csv(orig / "longitudinal.pairs.tsv", sep="\t"))
    n_l = unordered(pd.read_csv(next(nf.glob("contamination/*.longitudinal.pairs.tsv")), sep="\t"))
    errors += cmp_frame("longitudinal", o_l, n_l, ["sample_a", "sample_b"],
                        ["hom_discordance", "n_both_hom", "somalier_relatedness", "somalier_concordance", "flag"])
    o_f = unordered(pd.read_csv(orig / "relatedness_flags.tsv", sep="\t"))
    n_f = unordered(pd.read_csv(next(nf.glob("contamination/*.relatedness_flags.tsv")), sep="\t"))
    errors += cmp_frame("relatedness_flags", o_f, n_f, ["sample_a", "sample_b"],
                        ["flags", "hom_discordance", "somalier_relatedness", "somalier_concordance"])
    print("ALL EQUIVALENT" if not errors else f"{errors} DIFFERENCES")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
