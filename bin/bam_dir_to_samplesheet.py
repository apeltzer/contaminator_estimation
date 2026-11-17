#!/usr/bin/env python3
"""
Build an nf-core/contaminationqc samplesheet (sample,bam,bai,individual,batch) from BAM/CRAM files
in S3 prefixes and/or local directories.

Each location is listed once (`aws s3api list-objects-v2`, paginated; no per-file requests), and
every BAM/CRAM is paired with its index from that listing (x.bam.bai, x.bai, x.bam.csi, x.cram.crai,
x.crai). Files without an index are skipped with a warning.

Examples
  # all BAMs below a prefix, batch = name of the BAM's parent folder
  bam_dir_to_samplesheet.py s3://bucket/runs/run42/ --batch-level 1 -o samplesheet.csv

  # several runs, batch from a path component, patient ID from the sample name, skip intermediates
  bam_dir_to_samplesheet.py s3://bucket/runs/run42/ s3://bucket/runs/run43/ \\
      --batch-regex 'runs/([^/]+)/' --individual-regex '^(P\\d+)_' --exclude '_tmp|Logs_Intermediates' \\
      -o samplesheet.csv

Options for S3 access (--profile, --no-sign-request, --region) are passed to the AWS CLI.
Requires: python3, AWS CLI v2 (for s3:// locations).
"""
import argparse
import collections
import csv
import json
import os
import re
import subprocess
import sys

ALN = re.compile(r"\.(bam|cram)$")
INDEX_SUFFIXES = {"bam": [".bam.bai", ".bai", ".bam.csi", ".csi"], "cram": [".cram.crai", ".crai"]}
VALID = re.compile(r"[^A-Za-z0-9._-]")


def log(msg):
    print(f"[bam_dir_to_samplesheet] {msg}", file=sys.stderr)


def list_s3(url, aws_args):
    bucket, _, prefix = url[5:].partition("/")
    cmd = ["aws", "s3api", "list-objects-v2", "--bucket", bucket, "--prefix", prefix,
           "--query", "Contents[].Key", "--output", "json"] + aws_args
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"listing {url} failed:\n{r.stderr.strip()}")
    keys = json.loads(r.stdout or "null") or []
    return [f"s3://{bucket}/{k}" for k in keys]


def list_local(path):
    out = []
    for root, _dirs, files in os.walk(path, followlinks=True):
        out += [os.path.abspath(os.path.join(root, f)) for f in files]
    return out


def first_group(regex, text, what):
    if not regex:
        return ""
    m = re.search(regex, text)
    if not m:
        return ""
    if not m.groups():
        sys.exit(f"--{what}-regex needs a capture group: {regex}")
    return m.group(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("locations", nargs="+", help="s3://bucket/prefix/ and/or local directories")
    ap.add_argument("-o", "--output", help="output CSV (default: stdout)")
    ap.add_argument("--include", help="only keep alignment files whose path matches this regex")
    ap.add_argument("--exclude", help="drop alignment files whose path matches this regex")
    ap.add_argument("--sample-regex", help="sample ID = first capture group on the file name (default: name without extension)")
    ap.add_argument("--individual-regex", help="individual = first capture group on the sample ID")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--batch", help="fixed batch for all samples")
    g.add_argument("--batch-level", type=int, help="batch = Nth parent folder name of the file (1 = direct parent)")
    g.add_argument("--batch-regex", help="batch = first capture group on the full path")
    ap.add_argument("--prefix-parent", action="store_true",
                    help="make duplicated sample IDs unique by prefixing the parent folder name")
    ap.add_argument("--profile", help="AWS CLI profile")
    ap.add_argument("--region", help="AWS region")
    ap.add_argument("--no-sign-request", action="store_true", help="anonymous S3 access (public buckets)")
    args = ap.parse_args()

    aws_args = (["--profile", args.profile] if args.profile else []) + \
               (["--region", args.region] if args.region else []) + \
               (["--no-sign-request"] if args.no_sign_request else [])

    files = []
    for loc in args.locations:
        listed = list_s3(loc, aws_args) if loc.startswith("s3://") else list_local(loc)
        log(f"{loc}: {len(listed)} objects")
        files += listed
    present = set(files)

    rows, skipped = [], 0
    for f in sorted(files):
        m = ALN.search(f)
        if not m:
            continue
        if args.include and not re.search(args.include, f):
            continue
        if args.exclude and re.search(args.exclude, f):
            continue
        stem = f[: m.start()]
        idx = next((stem + s for s in INDEX_SUFFIXES[m.group(1)] if stem + s in present), None)
        if not idx:
            log(f"WARNING: no index for {f} - skipped")
            skipped += 1
            continue
        parts = f.rstrip("/").split("/")
        name = ALN.sub("", parts[-1])
        sample = first_group(args.sample_regex, parts[-1], "sample") if args.sample_regex else name
        if args.batch is not None:
            batch = args.batch
        elif args.batch_level:
            batch = parts[-1 - args.batch_level] if len(parts) > args.batch_level else ""
        else:
            batch = first_group(args.batch_regex, f, "batch")
        rows.append({"sample": VALID.sub("_", sample), "bam": f, "bai": idx,
                     "individual": "", "batch": VALID.sub("_", batch), "_parent": VALID.sub("_", parts[-2])})

    counts = collections.Counter(r["sample"] for r in rows)
    dups = sorted(s for s, n in counts.items() if n > 1)
    if dups and args.prefix_parent:
        for r in rows:
            if r["sample"] in dups:
                r["sample"] = f"{r['_parent']}_{r['sample']}"
        counts = collections.Counter(r["sample"] for r in rows)
        dups = sorted(s for s, n in counts.items() if n > 1)
    if dups:
        sys.exit(f"duplicated sample IDs (use --prefix-parent, --sample-regex or --exclude): {', '.join(dups[:20])}"
                 + (" ..." if len(dups) > 20 else ""))
    for r in rows:
        r["individual"] = VALID.sub("_", first_group(args.individual_regex, r["sample"], "individual"))
    if not rows:
        sys.exit("no BAM/CRAM files with index found")

    out = open(args.output, "w", newline="") if args.output else sys.stdout
    w = csv.DictWriter(out, fieldnames=["sample", "bam", "bai", "individual", "batch"], extrasaction="ignore",
                       lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    if args.output:
        out.close()
    log(f"{len(rows)} samples written{f' to {args.output}' if args.output else ''}"
        f"{f', {skipped} skipped without index' if skipped else ''}; "
        f"{len({r['batch'] for r in rows})} batch(es), {sum(1 for r in rows if r['individual'])} with individual")


if __name__ == "__main__":
    main()
