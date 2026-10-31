#!/usr/bin/env bash
# mappability.sh -- genome in, mappability track out.  Works for any read type.
#
#   bash mappability.sh -g genome.fa -o outprefix [options]
#
#   -g FILE    genome FASTA (required)
#   -o PREFIX  output prefix (required); writes PREFIX.mappability.bedgraph, PREFIX.params.json, PREFIX.log
#   -B FILE    learn read/fragment lengths from this BAM (ATAC, WGS, HiFi, ONT ...). Without it:
#              Bubb et al. 2025 ATAC-style default (50 bp paired, log-normal fragments, median 200 bp)
#   -a NAME    aligner: auto (default) | bwa | minimap2
#   -p PRESET  minimap2 preset for long reads: map-hifi (default) | map-ont | map-pb
#   -s N       stride between fragment starts in bp (default: 10 for short reads, 500 for long reads)
#   -b N       window size in bp (default 100)
#   -c REGEX   chromosomes to simulate (default '^chr(0?[1-9]|10)$'; use '.' for all)
#   -t N       threads (default $SLURM_CPUS_PER_TASK or 8)
#   any other arguments are passed to simreads.py (e.g. --mode single --readlen 150 --sample 1000000)
#
# Steps: 1. simulate exact reads shaped like the data (simreads.py)
#        2. map them back (bwa mem for short reads, minimap2 for long reads)
#        3. per-window fraction correctly mapped with MAPQ > 0 (score_mappability.py)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THREADS=${SLURM_CPUS_PER_TASK:-8}; STRIDE=""; BIN=100; CHROMS='^chr(0?[1-9]|10)$'
GENOME=""; OUT=""; FROMBAM=""; ALIGNER="auto"; PRESET="map-hifi"; EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g) GENOME="$2"; shift 2;;   -o) OUT="$2"; shift 2;;     -B) FROMBAM="$2"; shift 2;;
    -a) ALIGNER="$2"; shift 2;;  -p) PRESET="$2"; shift 2;;  -s) STRIDE="$2"; shift 2;;
    -b) BIN="$2"; shift 2;;      -c) CHROMS="$2"; shift 2;;  -t) THREADS="$2"; shift 2;;
    *) EXTRA+=("$1"); shift;;
  esac
done
[[ -n "$GENOME" && -n "$OUT" ]] || { sed -n '2,20p' "$0"; exit 1; }
mkdir -p "$(dirname "$OUT")"; LOG="$OUT.log"; : > "$LOG"; PAR="$OUT.params.json"

# ---- decide read shape first (cheap: samples the BAM), so we know the aligner and stride ----
SIMARGS=(--chroms "$CHROMS" --write-params "$PAR" "${EXTRA[@]}")
[[ -n "$FROMBAM" ]] && SIMARGS+=(--from-bam "$FROMBAM")
python3 "$HERE/simreads.py" "$GENOME" "${SIMARGS[@]}" --stride 100000000 --chroms '^$' >/dev/null 2>>"$LOG" || true
MODE=$(python3 -c "import json;print(json.load(open('$PAR'))['mode'])")
MEDLEN=$(python3 -c "import json;print(json.load(open('$PAR'))['median_read_length'])")
if [[ "$ALIGNER" == "auto" ]]; then ALIGNER=$([[ "$MEDLEN" -ge 500 ]] && echo minimap2 || echo bwa); fi
if [[ -z "$STRIDE" ]]; then STRIDE=$([[ "$MEDLEN" -ge 500 ]] && echo 500 || echo 10); fi
echo "[mappability] genome=$GENOME from_bam=${FROMBAM:-none} mode=$MODE median_read_len=$MEDLEN aligner=$ALIGNER stride=$STRIDE bin=$BIN threads=$THREADS chroms=$CHROMS" | tee -a "$LOG"

# ---- index ----
if [[ "$ALIGNER" == "bwa" ]]; then
  command -v bwa >/dev/null || { echo "bwa not on PATH"; exit 1; }
  [[ -e "$GENOME.bwt" ]] || { echo "[mappability] bwa index (one time, ~1 h for 2 Gb)" | tee -a "$LOG"; bwa index "$GENOME" 2>>"$LOG"; }
  if [[ "$MODE" == "paired" ]]; then MAP=(bwa mem -p -t "$THREADS" -v 1 "$GENOME" -); else MAP=(bwa mem -t "$THREADS" -v 1 "$GENOME" -); fi
else
  command -v minimap2 >/dev/null || { echo "minimap2 not on PATH"; exit 1; }
  MMI="${GENOME%.*}.$PRESET.mmi"
  [[ -e "$MMI" ]] || { echo "[mappability] minimap2 index (one time, ~10 min)" | tee -a "$LOG"; minimap2 -x "$PRESET" -d "$MMI" "$GENOME" 2>>"$LOG"; }
  MAP=(minimap2 -ax "$PRESET" -t "$THREADS" --secondary=no "$MMI" -)
fi

# ---- simulate | map | score ----
echo "[mappability] simulate -> ${MAP[0]} -> score" | tee -a "$LOG"
python3 "$HERE/simreads.py" "$GENOME" "${SIMARGS[@]}" --stride "$STRIDE" 2>>"$LOG" \
  | "${MAP[@]}" 2>>"$LOG" \
  | python3 "$HERE/score_mappability.py" --bin "$BIN" --params "$PAR" -o "$OUT.mappability.bedgraph" 2>>"$LOG"

grep -h "^\[simreads\] learned\|^\[simreads\] manual\|^\[score\]" "$LOG" | tail -n 3
awk '$6!="NA"{n++; s+=$6; if($6<0.8) low++} END{printf "[mappability] windows %d | mean mappability %.3f | <80%%: %.1f%%\n", n, s/n, 100*low/n}' "$OUT.mappability.bedgraph" | tee -a "$LOG"
echo "[mappability] done -> $OUT.mappability.bedgraph"
