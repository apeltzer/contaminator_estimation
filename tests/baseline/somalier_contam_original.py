#!/usr/bin/env python3
"""
somalier_contam.py - relatedness, contamination and contamination-SOURCE attribution for
targeted-panel BAMs (written with TSO500 v2 in mind) on top of somalier >= 0.3.2.

Steps
  1. somalier extract        per BAM/CRAM (parallel, cached; duplicate SM tags get auto-prefixed)
  2. somalier relate         relatedness, CNV-robust concordance, p_middling_ab, CHARR
  3. somalier contamination  CHARR per sample + anchored (conpair-like) all-vs-all MLEs.
                             NB: somalier's pair model treats the ANCHOR as the receiver's own
                             clean genotype and marginalises the contaminant over population AF,
                             i.e. it answers "how much foreign DNA", not "whose".
  4. source attribution      (this script, reads the .somalier allele counts directly)
                             at sites where receiver R is homozygous-like, compare
                               H_src: R = (1-a)*R_self + a*C           (C = candidate's genotypes)
                               H_unk: R = (1-a)*R_self + a*random      (HWE on sites-file AF)
                             LLR = max_a LL_src - max_a LL_unk on the same sites.
                             LLR >> 0: C explains R's foreign alleles better than a random person.
  5. identity                hom-hom discordance between samples. Robust to contamination of either
                             sample, so usable for longitudinal same-patient checks and swaps.

Outputs (--outdir)
  summary.samples.tsv       per sample: depth, p_middling_ab, CHARR, source-agnostic alpha, best source, call
  attribution.pairs.tsv     receiver x candidate: alpha, LLRs, informative sites, identity + somalier metrics
  longitudinal.pairs.tsv    pairs expected to be the same individual (sample sheet) or genotype-identical
  relatedness_flags.tsv     identity mismatches / unexpected identity / elevated relatedness
  attribution_llr.png       LLR heatmap (if matplotlib is installed)
  somalier/                 raw relate (+ html) and contamination outputs

Example
  python somalier_contam.py --fasta hg19.fa --sites sites.hg19.vcf.gz \\
      --panel-bed tso500v2_targets.bed --sample-sheet samples.tsv -o contam_out -j 8
  samples.tsv: tab-separated with header; column 'bam' required, 'individual' optional (patient ID).

Requires: somalier >= 0.3.2 on PATH (or --somalier), python3 with numpy + pandas.
"""
import argparse
import bisect
import collections
import concurrent.futures as cf
import gzip
import hashlib
import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

XY = {"X", "chrX", "NC_000023.10", "NC_000023.11", "Y", "chrY", "NC_000024.9", "NC_000024.10"}
# alpha grid: fine below 5 % where it matters, coarser above
GRID = np.unique(np.r_[np.arange(0, 0.05, 0.0005), np.arange(0.05, 0.5 + 1e-9, 0.005)])


def log(msg):
    print(f"[somalier_contam] {msg}", file=sys.stderr, flush=True)


def openx(path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def run(cmd, env=None):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, env=env)
    if r.returncode != 0:
        sys.exit(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}\n{r.stderr[-3000:]}")
    return r


# ----------------------------------------------------------------------------- inputs
def read_inputs(args):
    rows = []
    if args.sample_sheet:
        sh = pd.read_csv(args.sample_sheet, sep="\t", dtype=str).fillna("")
        if "bam" not in sh.columns:
            sys.exit("sample sheet needs a 'bam' column")
        for _, r in sh.iterrows():
            rows.append({"bam": r["bam"], "individual": r.get("individual", "")})
    rows += [{"bam": b, "individual": ""} for b in args.bams]
    if len(rows) < 2:
        sys.exit("need at least 2 BAM/CRAM files")
    for r in rows:
        if not os.path.exists(r["bam"]):
            sys.exit(f"missing input: {r['bam']}")
    return rows


def bam_contigs(path):
    """Reference names from a BAM header (pure python, BGZF is gzip-compatible). None for CRAM/errors."""
    if not str(path).endswith(".bam"):
        return None
    try:
        with gzip.open(path, "rb") as fh:
            if fh.read(4) != b"BAM\x01":
                return None
            (l_text,) = struct.unpack("<i", fh.read(4))
            fh.read(l_text)
            (n_ref,) = struct.unpack("<i", fh.read(4))
            names = set()
            for _ in range(n_ref):
                (l_name,) = struct.unpack("<i", fh.read(4))
                names.add(fh.read(l_name)[:-1].decode())
                fh.read(4)
            return names
    except Exception:
        return None


def site_contigs(path):
    seen = set()
    with openx(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                seen.add(line.split("\t", 1)[0])
    return seen


def restrict_sites(sites, bed, pad, out):
    """Keep only sites inside (padded) panel targets. Writes a plain-text VCF somalier can read."""
    if out.exists() and out.stat().st_mtime > max(os.path.getmtime(sites), os.path.getmtime(bed)):
        return out
    iv = collections.defaultdict(list)
    with openx(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split("\t")[:3]
            iv[c].append((max(0, int(s) - pad), int(e) + pad))
    merged = {}
    for c, lst in iv.items():
        lst.sort()
        m = []
        for s, e in lst:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[c] = ([x[0] for x in m], [x[1] for x in m])
    n_in = n_kept = 0
    tmp = out.with_suffix(".tmp")
    with openx(sites) as fh, open(tmp, "w") as o:
        for line in fh:
            if line.startswith("#"):
                o.write(line)
                continue
            n_in += 1
            c, p = line.split("\t", 2)[:2]
            if c in merged:
                st, en = merged[c]
                i = bisect.bisect_right(st, int(p) - 1) - 1
                if i >= 0 and int(p) - 1 < en[i]:
                    o.write(line)
                    n_kept += 1
    if n_kept == 0:
        sys.exit("no sites fall inside the panel BED - check contig naming (chr1 vs 1) and genome build")
    tmp.rename(out)
    log(f"panel restriction: kept {n_kept}/{n_in} sites")
    if n_kept < 300:
        log("WARNING: <300 sites in panel; consider `somalier find-sites --include panel.bed --snp-dist 2000`")
    return out


def read_sites_af(path):
    """Autosomal sites in somalier's storage order with AF of the stored B allele."""
    recs = []
    with openx(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if t[0] in XY:
                continue
            af = next((float(kv[3:].split(",")[0]) for kv in t[7].split(";") if kv.startswith("AF=")), None)
            if af is None:
                sys.exit(f"sites file lacks INFO/AF (needed for contamination): {path}")
            if not t[3] < t[4]:  # somalier stores alleles alphabetically -> flip AF like stored_alt_af()
                af = 1.0 - af
            recs.append((t[0], int(t[1]), af))
    recs.sort(key=lambda r: (r[0], r[1]))
    return np.clip(np.array([r[2] for r in recs], dtype=float), 1e-3, 1 - 1e-3)


# ----------------------------------------------------------------------------- somalier
def read_somalier(path):
    b = Path(path).read_bytes()
    if b[0] != 2:
        sys.exit(f"{path}: .somalier format version {b[0]}, expected 2 (somalier >= 0.3)")
    n = b[1]
    name = b[2:2 + n].decode()
    o = 2 + n
    na, _nx, _ny = struct.unpack_from("<HHH", b, o)
    o += 6
    c = np.frombuffer(b, dtype="<u4", count=3 * na, offset=o).reshape(na, 3).astype(float)
    return name, c  # columns: A-allele, B-allele, other


def extract_one(somalier, bam, sites, fasta, d, prefix=None):
    d.mkdir(parents=True, exist_ok=True)
    have = list(d.glob("*.somalier"))
    newest_input = max(os.path.getmtime(bam), os.path.getmtime(sites))
    if len(have) == 1 and have[0].stat().st_mtime >= newest_input:
        if prefix is None or read_somalier(have[0])[0].startswith(prefix):
            return have[0]
    for f in have:
        f.unlink()
    cmd = [somalier, "extract", "-d", d, "-s", sites, "-f", fasta]
    if prefix:
        cmd += ["--sample-prefix", prefix]
    run(cmd + [bam])
    out = list(d.glob("*.somalier"))
    if len(out) != 1:
        sys.exit(f"expected exactly one .somalier for {bam}, found {out}")
    return out[0]


def extract_all(args, rows, sites, outdir):
    def ddir(bam):
        h = hashlib.md5(os.path.abspath(bam).encode()).hexdigest()[:8]
        return outdir / "extract" / f"{Path(bam).stem}_{h}"

    with cf.ThreadPoolExecutor(args.threads) as ex:
        files = list(ex.map(lambda r: extract_one(args.somalier, r["bam"], sites, args.fasta, ddir(r["bam"])), rows))
    names = [read_somalier(f)[0] for f in files]
    dups = {n for n, k in collections.Counter(names).items() if k > 1}
    if dups:
        log(f"duplicate SM names {sorted(dups)} -> re-extracting with BAM-stem prefixes")
        stems = collections.Counter(Path(r["bam"]).stem for r in rows)
        for i, r in enumerate(rows):
            if names[i] in dups:
                st = Path(r["bam"]).stem
                pre = f"{st}-" if stems[st] == 1 else f"{st}_{i}-"
                files[i] = extract_one(args.somalier, r["bam"], sites, args.fasta, ddir(r["bam"]), pre)
                names[i] = read_somalier(files[i])[0]
    for r, f, n in zip(rows, files, names):
        r["somalier_file"], r["sample"] = str(f), n
    return rows


def run_somalier_cohort(args, rows, sites, outdir):
    sd = outdir / "somalier"
    sd.mkdir(exist_ok=True)
    files = [r["somalier_file"] for r in rows]
    env = dict(os.environ, SOMALIER_REPORT_ALL_PAIRS="1", SOMALIER_SEED="42")
    cmd = [args.somalier, "relate", "-s", sites, "-d", args.min_depth, "-o", sd / "relate"]
    groups = collections.defaultdict(list)
    for r in rows:
        if r["individual"]:
            groups[r["individual"]].append(r["sample"])
    groups = [g for g in groups.values() if len(g) > 1]
    if groups:
        (sd / "groups.txt").write_text("".join(",".join(g) + "\n" for g in groups))
        cmd += ["-g", sd / "groups.txt"]
    log("somalier relate")
    run(cmd + files, env=env)
    log("somalier contamination (all-vs-all)")
    run([args.somalier, "contamination", "-s", sites, "--min-depth", args.min_depth,
         "-o", sd / "contam"] + files)

    def rd(p):
        df = pd.read_csv(p, sep="\t")
        return df.rename(columns=lambda c: c.lstrip("#"))

    return (rd(sd / "relate.samples.tsv"), rd(sd / "relate.pairs.tsv"),
            rd(sd / "contam.samples.tsv"), rd(sd / "contam.pairs.tsv"))


# ----------------------------------------------------------------------------- models
def logsumexp0(x):
    m = x.max(0)
    return m + np.log(np.exp(x - m).sum(0))


def genotype_calls(C, min_depth, hom_max_af, het_min_af):
    ref, alt, oth = C[..., 0], C[..., 1], C[..., 2]
    dp = ref + alt
    ok = (dp >= min_depth) & (oth <= 0.1 * (dp + oth))
    af = np.divide(alt, dp, out=np.full_like(dp, np.nan), where=dp > 0)
    g = np.full(dp.shape, -1, dtype=np.int8)
    g[ok & (af <= hom_max_af)] = 0
    g[ok & (af >= 1 - hom_max_af)] = 2
    g[ok & (af >= het_min_af) & (af <= 1 - het_min_af)] = 1
    return g, alt, dp


def identity_matrix(G):
    """Hom-hom discordance: among sites where both samples are homozygous-like, the fraction with
    opposite homozygotes. ~0 for the same person, even if either sample carries a few % contamination."""
    H, A = (G == 0).astype(float), (G == 2).astype(float)
    both = (H + A) @ (H + A).T
    opp = H @ A.T + A @ H.T
    with np.errstate(invalid="ignore", divide="ignore"):
        return opp / both, both


def identity_groups(disc, nboth, max_disc, min_sites):
    n = len(disc)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if nboth[i, j] >= min_sites and disc[i, j] <= max_disc:
                parent[find(i)] = find(j)
    roots = {}
    return [roots.setdefault(find(i), f"geno{len(roots) + 1}") for i in range(n)]


def attribution(names, G, alt, dp, afB, args):
    N = len(names)
    eps = args.error_rate
    q = afB
    lpi = np.log(np.stack([(1 - q) ** 2, 2 * q * (1 - q), q ** 2]))  # 3 x S
    scale = np.minimum(1.0, args.max_eff_depth / np.maximum(dp, 1))  # temper over-confidence at TSO500 depth
    K, D = alt * scale, dp * scale
    a = GRID[:, None]
    pair_rows, self_rows = [], []
    for r in range(N):
        idx = np.where((G[r] == 0) | (G[r] == 2))[0]
        if len(idx) == 0:
            self_rows.append({"sample": names[r], "n_hom_sites": 0})
            continue
        h = G[r, idx][None, :] / 2.0
        k, n = K[r, idx][None, :], D[r, idx][None, :]
        M = np.empty((3, len(GRID), len(idx)))
        for g in (0, 1, 2):
            p = ((1 - a) * h + a * (g / 2.0)) * (1 - 2 * eps) + eps
            M[g] = k * np.log(p) + (n - k) * np.log1p(-p)
        Mu = logsumexp0(M + lpi[:, idx][:, None, :])  # contaminant genotype marginalised (HWE)
        llu = Mu.sum(1)
        self_rows.append({"sample": names[r], "n_hom_sites": len(idx),
                          "alpha_unknown_source": GRID[llu.argmax()],
                          "llr_contam_vs_clean": llu.max() - llu[0]})
        GC = G[:, idx].T  # S x N candidate genotypes at R's hom-like sites
        I = [(GC == g).astype(float) for g in (0, 1, 2)]
        LLs = M[0] @ I[0] + M[1] @ I[1] + M[2] @ I[2]  # grid x N
        LLu = Mu @ (I[0] + I[1] + I[2])
        iS, iU = LLs.argmax(0), LLu.argmax(0)
        cols = np.arange(N)
        lls, llu_c, ll0 = LLs[iS, cols], LLu[iU, cols], LLs[0]
        n_used = (GC >= 0).sum(0)
        n_inf = ((GC >= 0) & (GC != G[r, idx][:, None])).sum(0)
        for c in range(N):
            if c != r:
                pair_rows.append({"receiver": names[r], "candidate": names[c],
                                  "n_sites": int(n_used[c]), "n_informative": int(n_inf[c]),
                                  "alpha_source": GRID[iS[c]], "alpha_unknown_same_sites": GRID[iU[c]],
                                  "llr_vs_unknown": lls[c] - llu_c[c], "llr_vs_clean": lls[c] - ll0[c]})
    return pd.DataFrame(pair_rows), pd.DataFrame(self_rows)


# ----------------------------------------------------------------------------- report
def heatmap(pairs, names, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    m = pairs.pivot(index="receiver", columns="candidate", values="llr_vs_unknown").reindex(index=names, columns=names)
    pos = m.values[np.isfinite(m.values) & (m.values > 0)]
    lim = max(20.0, float(pos.max())) if pos.size else 20.0  # negatives saturate; sources stand out
    sz = max(5, 0.25 * len(names))
    fig, ax = plt.subplots(figsize=(sz + 2, sz))
    im = ax.imshow(m.values, cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax.set_xticks(range(len(names)), names, rotation=90, fontsize=6)
    ax.set_yticks(range(len(names)), names, fontsize=6)
    ax.set_xlabel("candidate source")
    ax.set_ylabel("receiver")
    fig.colorbar(im, label="LLR source vs unknown contaminant")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bams", nargs="*", help="BAM/CRAM files (alternatively/additionally --sample-sheet)")
    ap.add_argument("--sample-sheet", help="TSV with 'bam' and optional 'individual' (patient) columns")
    ap.add_argument("-f", "--fasta", required=True, help="reference FASTA the BAMs are aligned to")
    ap.add_argument("-s", "--sites", required=True, help="somalier sites VCF with INFO/AF, matching build/contig names")
    ap.add_argument("--panel-bed", help="restrict sites to panel targets (recommended for TSO500)")
    ap.add_argument("--bed-pad", type=int, default=100, help="pad BED intervals by this many bp [100]")
    ap.add_argument("-o", "--outdir", default="somalier_contam_out")
    ap.add_argument("--somalier", default="somalier", help="path to somalier >= 0.3.2 binary")
    ap.add_argument("-j", "--threads", type=int, default=4, help="parallel extract jobs [4]")
    ap.add_argument("--min-depth", type=int, default=20, help="min site depth (somalier recommends >=20 for cancer) [20]")
    ap.add_argument("--hom-max-af", type=float, default=0.15, help="minor-allele fraction tolerated for hom-like calls [0.15]")
    ap.add_argument("--het-min-af", type=float, default=0.25, help="het call window [x, 1-x] for candidate genotypes [0.25]")
    ap.add_argument("--error-rate", type=float, default=0.002, help="per-read error rate in the mixture likelihood [0.002]")
    ap.add_argument("--max-eff-depth", type=float, default=300, help="cap effective depth per site (overdispersion guard) [300]")
    ap.add_argument("--contam-threshold", type=float, default=0.01, help="flag samples with estimated contamination >= this [0.01]")
    ap.add_argument("--min-llr", type=float, default=10, help="LLR (nat) vs unknown contaminant to call a source [10]")
    ap.add_argument("--min-margin", type=float, default=5, help="LLR margin best vs 2nd-best individual [5]")
    ap.add_argument("--identity-max-disc", type=float, default=0.05, help="hom-hom discordance <= this = same individual [0.05]")
    ap.add_argument("--min-identity-sites", type=int, default=20, help="min shared hom sites for identity calls [20]")
    args = ap.parse_args()

    if not shutil.which(args.somalier):
        sys.exit(f"somalier not found: {args.somalier}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = read_inputs(args)

    # contig-naming sanity check (hg19 'chr1' vs GRCh37 '1' is the classic TSO500 trap)
    sc = site_contigs(args.sites)
    for r in rows:
        bc = bam_contigs(r["bam"])
        if bc is not None and not (bc & sc):
            sys.exit(f"no shared contig names between {r['bam']} and {args.sites} - wrong sites build/naming?")

    sites = args.sites
    if args.panel_bed:
        sites = str(restrict_sites(args.sites, args.panel_bed, args.bed_pad, outdir / "sites.panel.vcf"))

    log(f"somalier extract on {len(rows)} files ({args.threads} parallel)")
    rows = extract_all(args, rows, sites, outdir)
    rel_s, rel_p, con_s, con_p = run_somalier_cohort(args, rows, sites, outdir)

    # ---- load counts, genotype calls
    afB = read_sites_af(sites)
    names, counts = [], []
    for r in rows:
        n, c = read_somalier(r["somalier_file"])
        if len(c) != len(afB):
            sys.exit(f"{r['somalier_file']}: {len(c)} sites but sites file has {len(afB)} autosomal sites")
        names.append(n)
        counts.append(c)
    C = np.stack(counts)
    G, alt, dp = genotype_calls(C, args.min_depth, args.hom_max_af, args.het_min_af)
    log(f"{len(afB)} autosomal sites; median usable per sample: {int(np.median((G >= 0).sum(1)))}")

    disc, nboth = identity_matrix(G)
    geno = identity_groups(disc, nboth, args.identity_max_disc, args.min_identity_sites)
    ind = {r["sample"]: r["individual"] for r in rows}
    grp = dict(zip(names, geno))
    ix = {n: i for i, n in enumerate(names)}

    log("source attribution")
    pairs, selfdf = attribution(names, G, alt, dp, afB, args)

    # ---- annotate pairs with identity + somalier metrics
    pairs["hom_discordance"] = [disc[ix[a], ix[b]] for a, b in zip(pairs.receiver, pairs.candidate)]
    pairs["n_both_hom"] = [int(nboth[ix[a], ix[b]]) for a, b in zip(pairs.receiver, pairs.candidate)]
    pairs["same_genotype_group"] = [grp[a] == grp[b] for a, b in zip(pairs.receiver, pairs.candidate)]
    pairs["candidate_individual"] = pairs.candidate.map(ind)
    rp = rel_p.copy()
    rp["key"] = [frozenset(x) for x in zip(rp.sample_a, rp.sample_b)]
    rp = rp.set_index("key")[["relatedness", "concordance", "ibs0", "n"]].add_prefix("somalier_")
    pairs["key"] = [frozenset(x) for x in zip(pairs.receiver, pairs.candidate)]
    pairs = pairs.join(rp, on="key").drop(columns="key")
    cp = con_p.rename(columns={"sample_name": "receiver", "anchor_sample": "candidate",
                               "contamination_mle": "somalier_anchored_mle"})
    pairs = pairs.merge(cp[["receiver", "candidate", "somalier_anchored_mle"]], how="left")
    pairs = pairs.sort_values(["receiver", "llr_vs_unknown"], ascending=[True, False])
    pairs.to_csv(outdir / "attribution.pairs.tsv", sep="\t", index=False, float_format="%.4g")

    # ---- per-sample summary
    summ = []
    for r in rows:
        s = r["sample"]
        d = {"sample": s, "bam": r["bam"], "individual": r["individual"], "genotype_group": grp[s]}
        sub = pairs[(pairs.receiver == s) & ~pairs.same_genotype_group]
        if len(sub):
            best_by_grp = sub.assign(g=sub.candidate.map(grp)).sort_values("llr_vs_unknown", ascending=False)
            best_by_grp = best_by_grp.drop_duplicates("g")
            b = best_by_grp.iloc[0]
            second = best_by_grp.iloc[1].llr_vs_unknown if len(best_by_grp) > 1 else np.nan
            d.update(best_source=b.candidate, best_source_individual=b.candidate_individual,
                     source_alpha=b.alpha_source, source_llr_vs_unknown=b.llr_vs_unknown,
                     source_n_informative=b.n_informative, second_best_llr=second,
                     llr_margin=b.llr_vs_unknown - second if not np.isnan(second) else np.nan)
        same = pairs[(pairs.receiver == s) & pairs.same_genotype_group]
        d["somalier_anchored_mle_same_individual"] = same.somalier_anchored_mle.median() if len(same) else np.nan
        summ.append(d)
    summ = pd.DataFrame(summ).merge(selfdf, on="sample", how="left")
    keep = ["sample_id", "depth_mean", "gt_depth_mean", "n_het", "n_unknown", "p_middling_ab", "contamination_charr"]
    summ = summ.merge(rel_s[[c for c in keep if c in rel_s.columns]].rename(columns={"sample_id": "sample"}),
                      on="sample", how="left")

    def call(x):
        est = np.nanmax([x.get("contamination_charr", np.nan), x.get("alpha_unknown_source", np.nan)])
        llr, marg = x.get("source_llr_vs_unknown", np.nan), x.get("llr_margin", np.nan)
        if llr >= args.min_llr and (np.isnan(marg) or marg >= args.min_margin):
            c = "source_identified"
        elif llr >= args.min_llr:
            c = "source_ambiguous"
        elif est >= args.contam_threshold:
            c = "contaminated_source_not_in_cohort"
        else:
            c = "clean"
        if x.get("n_hom_sites", 0) < 50:
            c += ";low_sites"
        return c

    summ["call"] = summ.apply(call, axis=1)
    lead = ["sample", "individual", "genotype_group", "call", "contamination_charr", "alpha_unknown_source",
            "p_middling_ab", "best_source", "best_source_individual", "source_alpha", "source_llr_vs_unknown",
            "llr_margin", "source_n_informative", "somalier_anchored_mle_same_individual", "n_hom_sites"]
    summ = summ[[c for c in lead if c in summ.columns] + [c for c in summ.columns if c not in lead]]
    summ.to_csv(outdir / "summary.samples.tsv", sep="\t", index=False, float_format="%.4g")

    # ---- longitudinal / identity tables (unordered pairs)
    up = pairs[pairs.receiver < pairs.candidate].copy()
    up["individual_a"], up["individual_b"] = up.receiver.map(ind), up.candidate.map(ind)
    expected_same = (up.individual_a != "") & (up.individual_a == up.individual_b)
    known_diff = (up.individual_a != "") & (up.individual_b != "") & (up.individual_a != up.individual_b)
    long_ = up[expected_same | up.same_genotype_group].copy()
    long_["flag"] = np.where(expected_same[long_.index] & ~long_.same_genotype_group,
                             "IDENTITY_MISMATCH(swap/heavy contamination?)", "ok")
    cols = ["receiver", "candidate", "individual_a", "individual_b", "hom_discordance", "n_both_hom",
            "somalier_relatedness", "somalier_concordance", "flag"]
    long_[cols].rename(columns={"receiver": "sample_a", "candidate": "sample_b"}).to_csv(
        outdir / "longitudinal.pairs.tsv", sep="\t", index=False, float_format="%.4g")

    flags = []
    for _, x in up.iterrows():
        f = []
        if x.individual_a and x.individual_a == x.individual_b and not x.same_genotype_group:
            f.append("expected_same_but_discordant")
        if x.same_genotype_group and known_diff[_]:
            f.append("unexpected_identity(swap/duplicate)")
        if not x.same_genotype_group and x.somalier_relatedness > 0.2:
            f.append("elevated_relatedness(contamination/family)")
        if f:
            flags.append({"sample_a": x.receiver, "sample_b": x.candidate, "flags": ";".join(f),
                          "hom_discordance": x.hom_discordance, "somalier_relatedness": x.somalier_relatedness,
                          "somalier_concordance": x.somalier_concordance})
    pd.DataFrame(flags, columns=["sample_a", "sample_b", "flags", "hom_discordance", "somalier_relatedness",
                                 "somalier_concordance"]).to_csv(outdir / "relatedness_flags.tsv", sep="\t",
                                                                 index=False, float_format="%.4g")
    heatmap(pairs, names, outdir / "attribution_llr.png")

    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ[["sample", "call", "contamination_charr", "alpha_unknown_source", "best_source",
                    "source_alpha", "source_llr_vs_unknown", "llr_margin"]].to_string(index=False))
    log(f"done -> {outdir}/ (see somalier/relate.html for the interactive QC plots)")


if __name__ == "__main__":
    main()
