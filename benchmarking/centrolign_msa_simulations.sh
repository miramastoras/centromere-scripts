#!/bin/bash
#
# Slurm script to benchmark the output of a centrolign multiple sequence alignment
# of simulated sequences, created by the sim_centromere script.
# Calls the analyze_case.py and infer_tree.py scripts.
#
#
# Job name:
#SBATCH --job-name=jeizenga-simulated-centrolign-msa
#
# Partition - This is the queue it goes in:
#SBATCH --partition=medium
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
#SBATCH --mem=160gb
#
# Number of tasks (one for each CPU desired for use case):
#SBATCH --ntasks=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
#SBATCH --array=2-30
#SBATCH --output=array_job_%A_task_%a.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=12:00:00

date
hostname
pwd

CHR=X
DATE=20240308

CENTRODIR=/private/groups/patenlab/jeizenga/centromere
SIMDIR=$CENTRODIR/simulation/
CASEDIR=$SIMDIR/msa_chr"$CHR"_sim_cases_"$DATE"/

CASE=$(ls $CASEDIR | sort | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

source $CENTRODIR/venv/bin/activate
cd $CASEDIR/$CASE
mkdir -p induced
mkdir -p subprobs

CENTROLIGNDIR=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/
CENTROLIGN=$CENTROLIGNDIR/centrolign-3337cd7
TRUTH_COMPARE=$CENTROLIGNDIR/compare_truth_aln-3337cd7
TREE_COMPARE=$CENTROLIGNDIR/tree_compare-3337cd7
TREE_DIST=$CENTROLIGNDIR/tree_pair_dist-3337cd7
CENTRO_SCRIPTS_DIR=/private/groups/patenlab/jeizenga/GitHub/centromere-scripts/
ANALYZE_CASE=$CENTRO_SCRIPTS_DIR/analyze_case.py
INFER_TREE=$CENTRO_SCRIPTS_DIR/infer_tree.py

echo "beginning alignment of case" $CASEDIR/$CASE

/usr/bin/time -v $CENTROLIGN -v 3 -T tree.txt -A induced/aln -S subprobs/sp all_seqs.fasta > msa.gfa 2> >( tee err.txt >&2 )

echo "alignment completed, analyzing results"

/usr/bin/time -v $INFER_TREE induced 1 > inferred_tree.txt
/usr/bin/time -v $TREE_COMPARE tree.txt inferred_tree.txt > tree_comparison.tsv

# the output makes more sense if we give it a non-trivial directory
cd $CASEDIR
/usr/bin/time -v $ANALYZE_CASE $CASE $TREE_DIST $TRUTH_COMPARE

