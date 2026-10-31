# mappability

**Genome in, mappability track out, for any read type.**

A small, dependency-light pipeline that measures how confidently sequencing reads can be placed at
each position of a genome. It reproduces the short-read mappability method of
[Bubb et al. 2025, *Nature Plants*](https://doi.org/10.1038/s41477-025-02002-z) and extends it to
any read type: point it at a BAM you already have (ATAC-seq, Illumina WGS, PacBio HiFi / Fiber-seq,
ONT) and it simulates reads of that shape, maps them back, and reports the fraction placed correctly
per window.

Built in the [Ou Lab](https://github.com/oushujun) at The Ohio State University for filtering
Fiber-seq accessible chromatin regions (ACRs) in maize.

## Why

Maize is ~85% transposable elements. For 50 bp reads, half the genome is ambiguous; for 18 kb
Fiber-seq reads, almost none of it is. Which half of the genome you can trust depends entirely on
the reads you sequenced, so the mappability track should be built from reads like yours.

| B73 v5, chr1-10 | reads placed correctly | genome < 80% mappable | regular ACRs flagged (of 65,309) |
|---|---|---|---|
| 50 bp paired (paper style) | 68.0% | 49.0% | 3,967 |
| Fiber-seq-shaped, median 17.9 kb | 99.98% | 0.0% | 5 |

The only regions long reads cannot place are tandem arrays longer than a read: the rDNA cluster on
chr6 (16.79-17.22 Mb) and knob/satellite blocks on chr9 and chr7.

## Install

Requirements: Python 3.8+ with numpy; `bwa` (short reads); `minimap2` (long reads); `samtools`
(only for `-B`). On OSC Pitzer all four are in a conda env or module.

    git clone https://github.com/<you>/mappability.git
    cd mappability && bash test/run_toy_test.sh      # ~1 min, prints PASS

## Usage

    # paper-style: 50 bp paired reads, log-normal fragments (median 200 bp), bwa
    bash mappability.sh -g genome.fa -o out/genome_sr

    # shaped by your own data: learns read/fragment lengths, picks bwa or minimap2
    bash mappability.sh -g genome.fa -o out/genome_fiberseq -B my_fiberseq.bam
    bash mappability.sh -g genome.fa -o out/genome_atac     -B my_atac.bam

    # apply a track to peaks; flags regions below 80% (the paper's cutoff)
    python3 peak_mappability.py out/genome_sr.mappability.bedgraph peaks.bed -o peaks.mapp.bed

Output: `<prefix>.mappability.bedgraph` with columns `chrom start end correct total fraction`
(`NA` when no read was assigned to a window), plus `<prefix>.params.json` and `<prefix>.log`.

### Options

| flag | meaning | default |
|---|---|---|
| `-g` | genome FASTA | required |
| `-o` | output prefix | required |
| `-B` | BAM to learn read/fragment lengths from | none (paper default) |
| `-a` | aligner: `auto`, `bwa`, `minimap2` | `auto` (by median read length) |
| `-p` | minimap2 preset for long reads | `map-hifi` |
| `-s` | stride between fragment starts (bp); 1 = paper scale | 10 short / 500 long |
| `-b` | window size (bp) | 100 |
| `-c` | regex of chromosomes to simulate | `^chr(0?[1-9]\|10)$` |
| `-t` | threads | `$SLURM_CPUS_PER_TASK` or 8 |

Extra arguments pass through to `simreads.py`, e.g. `--mode single --readlen 150`,
`--frag-median 250 --frag-sigma 0.4`, `--sample 1000000`. Run `python3 simreads.py -h`.

## How it works

1. **Simulate** (`simreads.py`). Walk each chromosome; every `-s` bp cut a fragment with a length
   drawn from the learned (or default) distribution; take exact reads off it (both ends for paired,
   whole fragment for single, half reverse-complemented). No errors. True coordinates in the read
   name. Streamed as FASTQ.
2. **Map** (`bwa mem -p` or `minimap2 -ax <preset>`) back to the same genome.
3. **Score** (`score_mappability.py`). Correct = true chromosome, start within 5 bp, MAPQ > 0.
   Short reads count toward the window their outer end falls in (paper definition); long reads
   count toward every window they cover. Per window: `fraction = correct / total`.

See [docs/method.md](docs/method.md) for the paper's text and the design choices.

## Cluster examples

`slurm/run_short_read_B73.slurm` (paper-style, ~3 h on 40 cores incl. bwa index) and
`slurm/run_fiberseq_B73.slurm` (long-read, ~2.5 h) are the jobs used for the table above.
Indexes (`bwa index`, `minimap2 -d`) are built once next to the FASTA.

## Citation

If you use this, please cite Bubb et al. 2025 for the method and this repository for the
implementation. Devin Li, Ou Lab, The Ohio State University, 2026.

## License

MIT
