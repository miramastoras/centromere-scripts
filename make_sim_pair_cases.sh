#!/bin/bash
# Job name:
#SBATCH --job-name=jeizenga-make-pair-sim-cases
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
#SBATCH --mem=3gb
#
# Number of tasks (one for each CPU desired for use case):
#SBATCH --ntasks=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
#SBATCH --array=1-60
#SBATCH --output=array_job_%A_task_%a.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=1:00:00

date
hostname
pwd
echo "array job" $SLURM_ARRAY_TASK_ID

# note: sample size is handled in sbatch array size
CHR=X
DATE=20231215

CENTRO_DIR=/private/groups/patenlab/jeizenga/centromere
SIMDIR=$CENTRO_DIR/simulation/
OUTPARDIR=$SIMDIR/pair_chr"$CHR"_sim_cases_"$DATE"/

SIM_CENTROMERE=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/sim_centromere-429efbd

# the base array that we'll simulate
FASTA=$SIMDIR/chm13_chr"$CHR"_active_array.upper.fasta
BED=$SIMDIR/chr"$CHR"_shifted_hors.bed

GENS=("25" "50" "100" "150" "200", "300")
i=$(($SLURM_ARRAY_TASK_ID % 6))
GEN=${GENS[$i]}
OUTDIR=$OUTPARDIR/gen"$GEN"/case_"$SLURM_ARRAY_TASK_ID"/
mkdir -p $OUTDIR
# simulate the sequence
echo "generating sequence pair" $SLURM_ARRAY_TASK_ID "with" $GEN "generations"
/usr/bin/time -v $SIM_CENTROMERE -g $GEN -o $OUTDIR/sim $FASTA $BED
# record this case for the centrolign file
echo gen"$GEN"/case_"$SLURM_ARRAY_TASK_ID" >> $OUTPARDIR/cases.txt

# this is silly, but it makes the evaluation logic much simpler if we always have a tree
printf "(seq1,seq2):%d;\n" $GEN > $OUTDIR/tree.txt
cat $OUTDIR/sim*.fasta > $OUTDIR/all_seqs.fasta

