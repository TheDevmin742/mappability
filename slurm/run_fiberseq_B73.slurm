#!/bin/bash
#SBATCH --job-name=mapp_long
#SBATCH --account=PAS2387
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=120G
#SBATCH --output=mapp_long_%j.out
# Long-read (Fiber-seq-shaped) mappability of B73 v5: read lengths learned from PS00291,
# simulated every 500 bp, mapped with minimap2 map-hifi, every covered window scored.
set -euo pipefail
export PATH="$HOME/.conda/envs/quast_env/bin:$PATH"
cd ~/data/mappability
GENOME=~/data/fiber-seq/Zm-B73-REFERENCE-NAM-5.0.fa
BAM=~/data/fiber-seq/PS00291bc2062.fiberseq.chr_fix.bam
OUT=/fs/scratch/PAS2387/devin/mappability/B73v5_fiberseq
mkdir -p "$(dirname $OUT)"
bash mappability.sh -g $GENOME -o $OUT -B $BAM -t $SLURM_CPUS_PER_TASK -s 500 -b 100 --sample 300000
cp $OUT.mappability.bedgraph $OUT.params.json $OUT.log ~/data/mappability/
