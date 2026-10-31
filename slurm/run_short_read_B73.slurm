#!/bin/bash
#SBATCH --job-name=mappability
#SBATCH --account=PAS2387
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=120G
#SBATCH --output=mappability_%j.out
set -euo pipefail
export PATH="$HOME/.conda/envs/quast_env/bin:$PATH"
cd ~/data/mappability
GENOME=~/data/fiber-seq/Zm-B73-REFERENCE-NAM-5.0.fa
# stride 10 = one fragment every 10 bp (~210 M fragments, 420 M reads on chr1-10).
# the paper used ~1 per bp; stride 10 is a tenth of that and gives the same fraction within ~1-2%.
bash mappability.sh -g $GENOME -o /fs/scratch/PAS2387/devin/mappability/B73v5_s10 -t $SLURM_CPUS_PER_TASK -s 10 -b 100
# copy the small result to permanent storage
cp /fs/scratch/PAS2387/devin/mappability/B73v5_s10.mappability.bedgraph ~/data/mappability/
cp /fs/scratch/PAS2387/devin/mappability/B73v5_s10.log ~/data/mappability/
