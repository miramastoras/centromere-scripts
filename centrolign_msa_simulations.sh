#!/bin/bash
# Job name:
#SBATCH --job-name=jeizenga-simulated-centrolign-msa
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
#SBATCH --mem=160gb
#
# Number of tasks (one for each CPU desired for use case):
#SBATCH --ntasks=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
#SBATCH --array=1-1
#SBATCH --output=array_job_%A_task_%a.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=24:00:00

date
hostname
pwd

CENTRODIR=/private/groups/patenlab/jeizenga/centromere
SIMDIR=$CENTRODIR/simulation/
CASEDIR=$SIMDIR/msa_chrX_sim_cases_20231205/

CASE=$(ls $CASEDIR | sort | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

source $CENTRODIR/venv/bin/activate
cd $CASEDIR/$CASE
mkdir -p induced
mkdir -p subprobs

CENTROLIGNDIR=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/
CENTROLIGN=$CENTROLIGNDIR/centrolign-833c9ac
TRUTH_COMPARE=$CENTROLIGNDIR/compare_truth_aln
TREE_COMPARE=$CENTROLIGNDIR/tree_compare
TREE_DIST=$CENTROLIGNDIR/tree_pair_dist
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

