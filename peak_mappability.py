#!/usr/bin/env python3
"""
peak_mappability.py -- per-peak mappability from the bedGraph made by score_mappability.py.

For each BED region: fraction = sum(correct) / sum(total) over the bins it
overlaps. Adds columns: correct, total, fraction, and a flag (LOW if fraction
< --cutoff, else OK). Chromosome names are normalised (chr01 == chr1).

    python3 peak_mappability.py mappability.bedgraph peaks.bed -o peaks.mapp.bed --cutoff 0.80
"""
import argparse, sys, bisect
from collections import defaultdict

def norm(c):
    c = c.strip()
    if c.startswith("chr0") and len(c) == 5 and c[4].isdigit():
        return "chr" + c[4]
    return c

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bedgraph"); ap.add_argument("peaks")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--cutoff", type=float, default=0.80)
    a = ap.parse_args()

    starts, corr, tot, binsz = defaultdict(list), defaultdict(list), defaultdict(list), None
    with open(a.bedgraph) as fh:
        for line in fh:
            c, s, e, k, t, _ = line.rstrip("\n").split("\t")
            c = norm(c); s, e = int(s), int(e)
            if binsz is None: binsz = e - s
            starts[c].append(s); corr[c].append(int(k)); tot[c].append(int(t))

    n = nlow = nnodata = 0
    with open(a.peaks) as fh, open(a.out, "w") as o:
        for line in fh:
            if not line.strip() or line[0] in "#tb" and (line.startswith("#") or line.startswith("track") or line.startswith("browser")):
                continue
            f = line.rstrip("\n").split("\t")
            c = norm(f[0]); s, e = int(f[1]), int(f[2])
            if c not in starts:
                o.write("\t".join(f[:3] + f[3:] + ["0", "0", "NA", "NODATA"]) + "\n"); nnodata += 1; n += 1; continue
            i = bisect.bisect_right(starts[c], s) - 1
            i = max(i, 0)
            K = T = 0
            while i < len(starts[c]) and starts[c][i] < e:
                K += corr[c][i]; T += tot[c][i]; i += 1
            if T == 0:
                flag, fr = "NODATA", "NA"; nnodata += 1
            else:
                frv = K / T; fr = f"{frv:.4f}"
                flag = "LOW" if frv < a.cutoff else "OK"
                if flag == "LOW": nlow += 1
            o.write("\t".join(f + [str(K), str(T), fr, flag]) + "\n"); n += 1
    print(f"[peaks] {n:,} regions | LOW (<{a.cutoff}) {nlow:,} ({100*nlow/max(n,1):.1f}%) | no data {nnodata:,}", file=sys.stderr)

if __name__ == "__main__":
    main()
