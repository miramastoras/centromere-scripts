#!/bin/bash
#
# Slurm script to run the sim_centromere script to generate multiple sequence
# alignment problems
#
# Job name:
#SBATCH --job-name=jeizenga-make-sim-cases
#
# Partition - This is the queue it goes in:
#SBATCH --partition=main
#
# Where to send email
#SBATCH --mail-user=joeizeng@gmail.com
#
#SBATCH --mail-type=ALL
#
# Number of nodes you need per job:
#SBATCH --nodes=1
#
# Memory needed for the jobs.
#SBATCH --mem=16gb
#
# Number of tasks (one for each CPU desired for use case):
#SBATCH --ntasks=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
#SBATCH --array=1-30
#SBATCH --output=array_job_%A_task_%a.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=6:00:00

date
hostname
pwd
echo "array job" $SLURM_ARRAY_TASK_ID

# note: sample size is handled in sbatch array size
CHR=X
DATE=20231215

CENTRO_DIR=/private/groups/patenlab/jeizenga/centromere
SIMDIR=$CENTRO_DIR/simulation/
OUTDIR=$SIMDIR/msa_chr"$CHR"_sim_cases_"$DATE"/case_"$SLURM_ARRAY_TASK_ID"/

GEN_TREE=/private/groups/patenlab/jeizenga/GitHub/centromere-scripts/generate_tree.py
SIM_CENTROMERE=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/sim_centromere-429efbd

# the base array that we'll simulate
FASTA=$SIMDIR/chm13_chr"$CHR"_active_array.upper.fasta
BED=$SIMDIR/chr"$CHR"_shifted_hors.bed
# number of sequences in each case
N_SEQS=8
# expected number of generations to the root of the tree
N_GENS=200

mkdir -p $OUTDIR

# simulate the tree
source $CENTRO_DIR/venv/bin/activate
TREE=$OUTDIR/tree.txt
echo "making tree with" $N_SEQS "sequences and expected height" $N_GENS "generations"
/usr/bin/time -v $GEN_TREE $N_SEQS $N_GENS > $TREE
# simulate the sequence
echo "generating sequences according to tree"
/usr/bin/time -v $SIM_CENTROMERE -o $OUTDIR/sim -T $TREE $FASTA $BED
cat $OUTDIR/sim*.fasta > $OUTDIR/all_seqs.fasta
