#!/usr/bin/env python3
"""
score_mappability.py -- turn aligned simulated reads into a mappability track.

Reads SAM from stdin (bwa mem or minimap2 output on simreads.py reads). For every primary
alignment it recovers the read's TRUE position from the read name and asks whether the
aligner put it back in the right place with MAPQ > 0.

Two ways to assign a read to windows (--assign):
  outer   the read counts toward the ONE window its true outer end falls in
          (R1: fragment start; R2: fragment end). This is the Bubb et al. 2025 definition
          for short reads: fraction_mapped = correct reads with MAPQ > 0 / all reads whose
          outer end falls in the window.
  cover   the read counts toward EVERY window it overlaps. Natural for long reads, where a
          sparse simulation still gives each window many covering reads; mappability of a
          window = fraction of reads covering it that were placed correctly.
  auto    (default) outer if the median read length < 500 bp, else cover.

Output (bedGraph-like, 6 columns):  chrom  start  end  correct  total  fraction
"""
import argparse, sys, json
import numpy as np

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--bin", type=int, default=100)
    ap.add_argument("--tol", type=int, default=5, help="bp tolerance on the alignment start")
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--assign", choices=["auto", "outer", "cover"], default="auto")
    ap.add_argument("--params", help="JSON from simreads.py --write-params (used by --assign auto)")
    a = ap.parse_args()

    assign = a.assign
    if assign == "auto":
        med = 50
        if a.params:
            try: med = json.load(open(a.params)).get("median_read_length", 50)
            except Exception: pass
        assign = "cover" if med >= 500 else "outer"
    print(f"[score] assign={assign} bin={a.bin} tol={a.tol} min_mapq={a.min_mapq}", file=sys.stderr)

    lengths, correct, total = {}, {}, {}
    n = n_ok = n_unm = n_lowq = n_wrong = 0
    B = a.bin

    def arrays(c, need):
        if c not in total:
            nb = (lengths.get(c, need + 1) + B - 1) // B
            total[c] = np.zeros(max(nb, 1), np.int32); correct[c] = np.zeros(max(nb, 1), np.int32)
        if need // B >= len(total[c]):
            grow = need // B + 1 - len(total[c])
            total[c] = np.concatenate([total[c], np.zeros(grow, np.int32)])
            correct[c] = np.concatenate([correct[c], np.zeros(grow, np.int32)])

    for line in sys.stdin:
        if line[0] == "@":
            if line.startswith("@SQ"):
                f = dict(x.split(":", 1) for x in line.rstrip("\n").split("\t")[1:])
                lengths[f["SN"]] = int(f["LN"])
            continue
        f = line.split("\t", 10)
        flag = int(f[1])
        if flag & 0x900: continue                 # secondary / supplementary
        qn = f[0]
        parts = qn.rsplit("_", 3) if qn.endswith(("_+", "_-")) else qn.rsplit("_", 2)
        if len(parts) == 4:                       # single-end: chrom_start_end_strand
            tchrom, ts, te, _ = parts; ts, te = int(ts), int(te)
            exp_pos, lo, hi = ts, ts, te
        else:                                     # paired: chrom_start_end (+ /1 /2 handled by flags)
            tchrom, ts, te = parts; ts, te = int(ts), int(te)
            L = len(f[9])
            if flag & 0x80: exp_pos, lo, hi = te - L, te - L, te      # R2 sits at the fragment end
            else:           exp_pos, lo, hi = ts, ts, ts + L
        outer = (te - 1) if (flag & 0x80) else ts
        arrays(tchrom, hi)
        if assign == "outer":
            b = outer // B; total[tchrom][b] += 1
        else:
            b0, b1 = lo // B, (hi - 1) // B; total[tchrom][b0:b1 + 1] += 1
        n += 1
        if flag & 0x4: n_unm += 1; continue
        if int(f[4]) < a.min_mapq: n_lowq += 1; continue
        if f[2] != tchrom or abs(int(f[3]) - 1 - exp_pos) > a.tol: n_wrong += 1; continue
        if assign == "outer": correct[tchrom][b] += 1
        else: correct[tchrom][b0:b1 + 1] += 1
        n_ok += 1

    with open(a.out, "w") as o:
        for c in sorted(total, key=lambda x: (len(x), x)):
            t, k = total[c], correct[c]; glen = lengths.get(c, len(t) * B)
            for i in range(len(t)):
                s = i * B; e = min(s + B, glen)
                if e <= s: break
                o.write(f"{c}\t{s}\t{e}\t0\t0\tNA\n" if t[i] == 0 else f"{c}\t{s}\t{e}\t{k[i]}\t{t[i]}\t{k[i]/t[i]:.4f}\n")
    print(f"[score] reads {n:,} | correct {n_ok:,} ({100*n_ok/max(n,1):.2f}%) | unmapped {n_unm:,} | MAPQ<{a.min_mapq} {n_lowq:,} | wrong place {n_wrong:,}", file=sys.stderr)

if __name__ == "__main__":
    main()
