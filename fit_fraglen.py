#!/usr/bin/env python3
"""
fit_fraglen.py -- fit a log-normal to ATAC-seq fragment lengths (optional step 0).

Reads a coordinate-sorted, paired-end ATAC BAM and fits ln(TLEN) ~ Normal(mu, sigma).
Prints --frag-median and --frag-sigma values to pass to simreads.py / mappability.sh.

    python3 fit_fraglen.py atac.bam [--max 2000] [--sample 2000000]
Needs pysam (available in the smapshot venv).
"""
import argparse, math, sys
import numpy as np, pysam
ap = argparse.ArgumentParser(); ap.add_argument("bam"); ap.add_argument("--max", type=int, default=2000); ap.add_argument("--sample", type=int, default=2000000)
a = ap.parse_args()
tl = []
for r in pysam.AlignmentFile(a.bam):
    if r.is_proper_pair and r.is_read1 and not r.is_secondary and not r.is_supplementary and 0 < r.template_length <= a.max:
        tl.append(r.template_length)
        if len(tl) >= a.sample: break
tl = np.array(tl, float); ln = np.log(tl)
mu, sig = ln.mean(), ln.std()
print(f"n={len(tl):,}  median={np.median(tl):.0f}  mean={tl.mean():.0f}")
print(f"log-normal fit: mu={mu:.3f} sigma={sig:.3f}  ->  --frag-median {math.exp(mu):.0f} --frag-sigma {sig:.2f}")
