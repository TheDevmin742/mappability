# Changelog

## v2.0 (2026-10-18)
- `-B example.bam`: learn read and fragment lengths from any BAM; simulate single or paired reads.
- Automatic aligner choice: bwa for short reads, minimap2 for long reads.
- `cover` window assignment for long reads; `outer` (paper) for short reads.
- `mappability.sh` writes `<prefix>.params.json`; summary line with mean and % windows < 80%.
- Self-test `test/run_toy_test.sh`.

## v1.0 (2026-10-10)
- Paper-style pipeline: 50 bp paired exact reads, log-normal fragments, bwa, outer-end scoring.
- `peak_mappability.py` to apply a track to a peak BED.
