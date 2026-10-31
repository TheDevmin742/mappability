# Method

## Definition
Mappability of a genomic window = the fraction of simulated reads originating in (short reads) or
covering (long reads) that window that an aligner places back at the correct position with MAPQ > 0.

## Source
Bubb, K.L., Hamm, M.O., Tullius, T.W. et al. The regulatory potential of transposable elements in
maize. *Nature Plants* 11, 1181-1192 (2025). https://doi.org/10.1038/s41477-025-02002-z
Methods, "Short-read mappability analysis":

> A total of 2.1 billion fragments were generated evenly distributed across the B73 reference genome
> chromosomes 1-10. Fragment lengths were sampled from a log-normal distribution fit to one of our
> ATAC-seq data sets. For each simulated fragment, a paired-end read was generated with 50 base reads
> on either end of the fragment. The true start and end of the fragment was encoded in the read name.
> We did not simulate per-base errors in these reads; each read matches exactly the reference sequence
> from which it was generated. These reads were then mapped back to the genome using BWA. The
> 'fraction mapped' for a given region or window was calculated as the number of correctly mapped
> reads with mapq score >0, divided by the total number of simulated reads with the outer end (Tn5
> insertion) falling in the region.

and "Classification of ACRs":

> (3) ACRs with mappability <80% and both high FIRE accessibility score and FE length >=200 bp were
> included in the 'Low mappability' set

## What this implementation adds
1. **Any read type.** `-B example.bam` samples the user's BAM for read lengths (and insert sizes if
   paired) and simulates reads from those empirical distributions. Without `-B`, the paper's default
   (50 bp paired, log-normal fragments, median 200 bp) is used.
2. **Aligner by read length.** bwa mem for reads < 500 bp (paired reads interleaved with `-p`),
   minimap2 for longer reads (`map-hifi` default; `-p map-ont` / `map-pb`).
3. **Two window assignments.** Short reads: outer-end rule as in the paper. Long reads: a read
   counts toward every window it covers, so a sparse simulation (one read per 500 bp) still gives
   ~35 reads per 100 bp window.
4. **Streaming.** Reads go simulator -> aligner -> scorer through pipes; nothing hits the disk but
   the final track.

## Correctness check
A read is "correct" when it maps to its true chromosome with alignment start within `--tol` (5 bp)
of its true start and MAPQ >= `--min-mapq` (1). Secondary and supplementary alignments are ignored.
Reads containing N are not simulated.

## Known limits
- No sequencing errors are simulated (as in the paper). Real reads score marginally lower.
- The sampled BAM is read from the top; for a sorted BAM that is the start of chromosome 1, which is
  fine for length distributions but would not be for anything position-dependent.
- Fragment and read lengths are drawn independently in paired mode.
