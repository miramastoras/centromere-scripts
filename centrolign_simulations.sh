#!/bin/bash
# Job name:
#SBATCH --job-name=jeizenga-simulated-centrolign
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
#SBATCH --mem=56gb
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
#SBATCH --time=1:00:00

date
hostname
pwd

# list cases in $SIMDIR/cases.txt

SIMDIR=/private/groups/patenlab/jeizenga/centromere/simulation/pair_chrX_sim_cases_20231214
CASE=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$SIMDIR"/cases.txt)
CASEDIR=$SIMDIR/$CASE
WORKDIR=$SIMDIR/work

CENTROLIGN_OUTFILE=$CASEDIR/aln_centrolign.txt
UNIALIGNER_OUTFILE=$CASEDIR/aln_unialigner.txt
WFA_OUTFILE=$CASEDIR/aln_wfa.txt

CENTROLIGN_DIR=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/
SCRIPTS_DIR=/private/groups/patenlab/jeizenga/GitHub/centromere-scripts/

CENTROLIGN=$CENTROLIGN_DIR/centrolign
TRUTH_COMPARE=$CENTROLIGN_DIR//compare_truth_aln
UNIALIGNER=/private/groups/patenlab/jeizenga/GitHub/unialigner/tandem_aligner/build/bin/tandem_aligner
TO_RAW_SEQ=$SCRIPTS_DIR/fasta_to_raw_seq.py
ANALYZE_CASE=$SCRIPTS_DIR/analyze_pair_case.py
WFA=/private/groups/patenlab/jeizenga/GitHub/WFA2-lib/bin/align_benchmark

# scores: M,X,O1,E1,O2,E2
# heuristic: min distance, max distance from lead, reduction interval
WFA_PARAMS="-a gap-affine2p-wfa --affine2p-penalties -20,80,100,30,5000,1 --wfa-heuristic wfa-adaptive --wfa-heuristic-parameters 1000,20000,50 --wfa-memory ultralow --wfa-span global"

mkdir -p $WORKDIR
cd $WORKDIR

echo "simulation:" $CASEDIR

FASTA1=${CASEDIR}/sim_seq1.fasta
FASTA2=${CASEDIR}/sim_seq2.fasta

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2

echo "aligning with centrolign"
TEMP_FASTA=${WORKDIR}/sim_joined_"$SLURM_ARRAY_TASK_ID".fa
mkdir -p `dirname $TEMP_FASTA`
cat $FASTA1 $FASTA2 > $TEMP_FASTA
mkdir -p `dirname $CENTROLIGN_OUTFILE`
/usr/bin/time -v ${CENTROLIGN} -v 3 $TEMP_FASTA > $CENTROLIGN_OUTFILE
rm $TEMP_FASTA

echo "aligning with unaligner"
mkdir -p `dirname $UNIALIGNER_OUTFILE`
UNIALIGNER_TEMP_OUTDIR=$WORKDIR/tmp_out_"$SLURM_ARRAY_TASK_ID"
/usr/bin/time -v ${UNIALIGNER} --first $FASTA1 --second $FASTA2 -o $UNIALIGNER_TEMP_OUTDIR
mv $UNIALIGNER_TEMP_OUTDIR/cigar.txt $UNIALIGNER_OUTFILE
# delete the rest of the output
rm -r $UNIALIGNER_TEMP_OUTDIR

echo "aligning with WFA"
RAW_SEQ_TEMP=${WORKDIR}/tmp_raw_seq_"$SLURM_ARRAY_TASK_ID".txt
WFA_TEMP_OUT=${WORKDIR}/tmp_wfa_out_"$SLURM_ARRAY_TASK_ID".txt
$TO_RAW_SEQ $FASTA1 > $RAW_SEQ_TEMP
$TO_RAW_SEQ $FASTA2 >> $RAW_SEQ_TEMP
# limit memory to 32 gB and runtime to 30 min, but don't consider it a failure if we don't get it
true || timeout -v 30m ulimit -m 33554432 $WFA $WFA_PARAMS -i $RAW_SEQ_TEMP -o $WFA_TEMP_OUT
# remove the score from the output
if [ -f $WFA_TEMP_OUT ]; then
    cut -f 2 $WFA_TEMP_OUT > $WFA_OUTFILE
else
    touch $WFA_OUTFILE
fi
rm -f $RAW_SEQ_TEMP
rm -f $WFA_TEMP_OUT

# do this from outside the directory to get more sensible output
cd $SIMDIR
$ANALYZE_CASE $CASE $TRUTH_COMPARE

