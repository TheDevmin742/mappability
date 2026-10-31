#!/usr/bin/env python3
"""
simreads.py -- simulate exact, evenly spaced reads from a genome, shaped like real data.

Two ways to set the read shape:
  --from-bam example.bam   learn read lengths (and insert sizes, if paired) from any BAM:
                           ATAC-seq, Illumina WGS, PacBio HiFi, ONT ... The simulation then
                           draws read / fragment lengths from those empirical distributions.
  manual                   --mode paired --readlen 50 --frag-median 200 --frag-sigma 0.45
                           (the Bubb et al. 2025 ATAC-style default), or --mode single --readlen N.

Fragments start every --stride bp along each chromosome (evenly distributed, as in the paper).
Reads are exact copies of the reference (no errors) and the true coordinates are in the name:
    paired:  @<chrom>_<start>_<end>/1  and  /2     (interleaved R1,R2 -> `bwa mem -p`)
    single:  @<chrom>_<start>_<end>_<strand>       (half of the reads are reverse-complemented)
Writes FASTQ to stdout. Prints the learned/used length summary to stderr, and with
--write-params FILE also writes a small JSON that mappability.sh reads to pick the aligner.
"""
import argparse, sys, re, math, json, subprocess, itertools
import numpy as np

def revcomp(s, _t=str.maketrans("ACGTacgtNn", "TGCAtgcaNn")):
    return s.translate(_t)[::-1]

def read_fasta(path, keep):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and keep(name):
                    yield name, "".join(buf).upper()
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if name is not None and keep(name):
            yield name, "".join(buf).upper()

def learn_from_bam(bam, nsample, samtools="samtools"):
    """Return dict(mode, readlens, fraglens) sampled from the first nsample primary reads."""
    cmd = [samtools, "view", "-F", "0x904", bam]   # primary, mapped, not secondary/supplementary
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    rl, fl, paired = [], [], 0
    for line in itertools.islice(p.stdout, nsample):
        f = line.split("\t", 10)
        flag = int(f[1]); seq = f[9]
        rl.append(len(seq))
        if flag & 0x1:
            paired += 1
            if (flag & 0x2) and (flag & 0x40):            # proper pair, read1 only (count each fragment once)
                t = abs(int(f[8]))
                if t > 0: fl.append(t)
    p.stdout.close(); p.wait()
    if not rl:
        sys.exit(f"[simreads] ERROR: no primary mapped reads found in {bam} (is it a valid, non-empty BAM?)")
    mode = "paired" if paired >= 0.5 * max(len(rl), 1) and fl else "single"
    return {"mode": mode, "readlens": np.array(rl), "fraglens": np.array(fl) if fl else None, "n": len(rl)}

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome")
    ap.add_argument("--from-bam", help="learn read/fragment length distributions from this BAM")
    ap.add_argument("--sample", type=int, default=500000, help="reads to sample from --from-bam")
    ap.add_argument("--mode", choices=["paired", "single"], help="override: paired-end or single reads")
    ap.add_argument("--stride", type=int, default=10, help="bp between fragment starts (1 = paper scale). Use ~500-1000 for long reads.")
    ap.add_argument("--readlen", type=int, default=50, help="manual read length (ignored with --from-bam)")
    ap.add_argument("--frag-median", type=float, default=200.0, help="manual log-normal fragment median (paired)")
    ap.add_argument("--frag-sigma", type=float, default=0.45, help="manual log-normal sigma (paired)")
    ap.add_argument("--frag-min", type=int, default=50)
    ap.add_argument("--frag-max", type=int, default=2000, help="manual cap for paired fragments")
    ap.add_argument("--max-len", type=int, default=100000, help="cap on any simulated read/fragment length")
    ap.add_argument("--chroms", default=r"^chr(0?[1-9]|10)$", help="regex of chromosomes to simulate")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=200000)
    ap.add_argument("--write-params", help="write a JSON summary (mode, median read length ...) here")
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    pat = re.compile(a.chroms)

    # ---- read shape -------------------------------------------------------------------------
    if a.from_bam:
        L = learn_from_bam(a.from_bam, a.sample, a.samtools)
        mode = a.mode or L["mode"]
        readlens = L["readlens"]; fraglens = L["fraglens"]
        print(f"[simreads] learned from {a.from_bam}: {L['n']:,} reads, mode={mode}, "
              f"read length median {int(np.median(readlens))} (p5 {int(np.percentile(readlens,5))}, p95 {int(np.percentile(readlens,95))})"
              + (f", fragment median {int(np.median(fraglens))}" if fraglens is not None and mode == "paired" else ""), file=sys.stderr)
        def draw_read(n): return np.minimum(rng.choice(readlens, size=n), a.max_len)
        def draw_frag(n): return np.minimum(rng.choice(fraglens, size=n), a.max_len)
        med_read = int(np.median(readlens))
    else:
        mode = a.mode or "paired"
        mu = math.log(a.frag_median)
        def draw_read(n): return np.full(n, a.readlen, dtype=np.int64)
        def draw_frag(n):
            fl = rng.lognormal(mu, a.frag_sigma, size=n)
            return np.clip(np.rint(fl), a.frag_min, a.frag_max).astype(np.int64)
        med_read = a.readlen
        print(f"[simreads] manual: mode={mode}, read length {a.readlen}"
              + (f", log-normal fragments median {a.frag_median:.0f} sigma {a.frag_sigma}" if mode == "paired" else ""), file=sys.stderr)

    if a.write_params:
        json.dump({"mode": mode, "median_read_length": int(med_read), "stride": a.stride,
                   "from_bam": a.from_bam}, open(a.write_params, "w"))

    out = sys.stdout
    n_frag = n_skip = 0
    for chrom, seq in read_fasta(a.genome, lambda n: pat.match(n) is not None):
        glen = len(seq)
        starts = np.arange(0, max(glen - 50, 1), a.stride, dtype=np.int64)
        for i in range(0, len(starts), a.chunk):
            st = starts[i:i + a.chunk]
            lines = []; ap_ = lines.append
            if mode == "paired":
                fl = draw_frag(len(st)).astype(np.int64)
                r1l = draw_read(len(st)).astype(np.int64); r2l = draw_read(len(st)).astype(np.int64)
                fl = np.maximum(fl, np.maximum(r1l, r2l))
                en = np.minimum(st + fl, glen)
                for s, e, l1, l2 in zip(st.tolist(), en.tolist(), r1l.tolist(), r2l.tolist()):
                    if e - s < max(l1, l2): continue
                    r1 = seq[s:s + l1]; r2 = revcomp(seq[e - l2:e])
                    if "N" in r1 or "N" in r2: n_skip += 1; continue
                    nm = f"{chrom}_{s}_{e}"
                    ap_(f"@{nm}/1\n{r1}\n+\n{'I'*l1}\n@{nm}/2\n{r2}\n+\n{'I'*l2}\n"); n_frag += 1
            else:
                rl = draw_read(len(st)).astype(np.int64)
                strand = rng.random(len(st)) < 0.5
                en = np.minimum(st + rl, glen)
                for s, e, rev in zip(st.tolist(), en.tolist(), strand.tolist()):
                    if e - s < 30: continue
                    r = seq[s:e]
                    if "N" in r: n_skip += 1; continue
                    if rev: r = revcomp(r)
                    ap_(f"@{chrom}_{s}_{e}_{'-' if rev else '+'}\n{r}\n+\n{'I'*(e-s)}\n"); n_frag += 1
            out.write("".join(lines))
        print(f"[simreads] {chrom}: {glen:,} bp done, cumulative {n_frag:,} {'fragments' if mode=='paired' else 'reads'} (skipped {n_skip:,} with N)", file=sys.stderr)
    print(f"[simreads] total {n_frag:,} {'fragments' if mode=='paired' else 'reads'}", file=sys.stderr)

if __name__ == "__main__":
    main()
