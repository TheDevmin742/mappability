#!/usr/bin/env bash
# Self-test: builds a 480 kb toy genome with a 20 kb repeat copied 5 times, runs the pipeline
# three ways (paper-style short reads; short reads learned from a BAM; long reads learned from a
# BAM) and checks that a peak inside the repeat is flagged LOW while a unique peak is OK.
# Needs: python3 + numpy, bwa, minimap2, samtools. Runs in ~1 minute.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; TOP="$(dirname "$HERE")"; OUT="$HERE/out"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
python3 - <<'PY'
import numpy as np
rng=np.random.default_rng(7); b=np.array(list("ACGT"))
def rnd(n): return "".join(b[rng.integers(0,4,n)])
rep=rnd(20000)
chr1=rnd(100000)+rep+rnd(50000)+rep+rnd(50000)+rep+rnd(30000)+rep+rnd(20000)+rep+rnd(30000)
chr2=rnd(100000)
with open("toy.fa","w") as f:
    for n,s in [("chr1",chr1),("chr2",chr2)]:
        f.write(f">{n}\n")
        for i in range(0,len(s),80): f.write(s[i:i+80]+"\n")
PY
printf "chr1\t50000\t50500\tunique_peak\nchr1\t105000\t105400\trepeat_peak\n" > peaks.bed
echo "== 1/3 paper-style 50 bp paired reads (bwa)"
bash "$TOP/mappability.sh" -g toy.fa -o sr -t 4 -s 20 -c '^chr[12]$' | tail -2
echo "== make toy BAMs"
python3 "$TOP/simreads.py" toy.fa --mode paired --readlen 75 --frag-median 300 --frag-sigma 0.3 --stride 300 --chroms '^chr[12]$' 2>/dev/null \
  | bwa mem -p -t 4 -v 1 toy.fa - 2>/dev/null | samtools sort -o short.bam - && samtools index short.bam
python3 "$TOP/simreads.py" toy.fa --mode single --readlen 8000 --stride 2000 --chroms '^chr[12]$' 2>/dev/null \
  | minimap2 -ax map-hifi -t 4 toy.fa - 2>/dev/null | samtools sort -o long.bam - && samtools index long.bam
echo "== 2/3 shaped by a short-read BAM"
bash "$TOP/mappability.sh" -g toy.fa -o short -t 4 -s 20 -c '^chr[12]$' -B short.bam | tail -2
echo "== 3/3 shaped by a long-read BAM (minimap2)"
bash "$TOP/mappability.sh" -g toy.fa -o long -t 4 -s 200 -c '^chr[12]$' -B long.bam | tail -2
echo "== per-peak"
ok=1
for r in sr short long; do
  python3 "$TOP/peak_mappability.py" $r.mappability.bedgraph peaks.bed -o $r.peaks.bed 2>/dev/null
  printf "%-6s %s\n" "$r" "$(cut -f4,7,8 $r.peaks.bed | tr '\t' ':' | tr '\n' ' ')"
  grep -q "unique_peak.*OK" $r.peaks.bed || ok=0
done
grep -q "repeat_peak.*LOW" sr.peaks.bed || ok=0
[[ $ok == 1 ]] && echo "PASS" || { echo "FAIL"; exit 1; }
