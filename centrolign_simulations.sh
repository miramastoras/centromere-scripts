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
#SBATCH --array=1-20
#SBATCH --output=array_job_%A_task_%a.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=1:00:00

date
hostname
pwd

CENTROLIGN=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/centrolign-f3438dd
UNIALIGNER=/private/groups/patenlab/jeizenga/GitHub/unialigner/tandem_aligner/build/bin/tandem_aligner
SIMDIR=/private/groups/patenlab/jeizenga/centromere/simulation/
WORKDIR=$SIMDIR/work/
OUTDIR=$SIMDIR/alignments_20231201/

mkdir -p $OUTDIR
mkdir -p $WORKDIR
cd $WORKDIR

SIM=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$SIMDIR"/simulations.txt)
echo "simulation:" $SIMDIR/$SIM
echo "out:" $OUTDIR

FASTA1=${SIMDIR}/${SIM}_seq1.fasta
FASTA2=${SIMDIR}/${SIM}_seq2.fasta

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2

echo "aligning with centrolign"
TEMP_FASTA=${WORKDIR}/${SIM}_joined.fa
mkdir -p `dirname $TEMP_FASTA`
cat $FASTA1 $FASTA2 > $TEMP_FASTA
CENTROLIGN_OUTFILE=$OUTDIR/"$SIM"_aln_centrolign.txt
mkdir -p `dirname $CENTROLIGN_OUTFILE`
${CENTROLIGN} -v 3 --skip-calibration $TEMP_FASTA > $CENTROLIGN_OUTFILE
rm $TEMP_FASTA

echo "aligning with unaligner"
UNIALIGNER_OUTFILE=$OUTDIR/"$SIM"_aln_unialigner.txt
mkdir -p `dirname $UNIALIGNER_OUTFILE`
UNIALIGNER_TEMP_OUTDIR=$WORKDIR/tmp_out_"$SLURM_ARRAY_TASK_ID"
${UNIALIGNER} --first $FASTA1 --second $FASTA2 -o $UNIALIGNER_TEMP_OUTDIR
mv $UNIALIGNER_TEMP_OUTDIR/cigar.txt $UNIALIGNER_OUTFILE
rm -r $UNIALIGNER_TEMP_OUTDIR

#contents of simulations.txt
#gen100/sim_z120_
#gen100/sim_z220_
#gen100/sim_z320_
#gen100/sim_z420_
#gen100/sim_z520_
#gen100/sim_z620_
#gen100/sim_z720_
#gen100/sim_z820_
#gen100/sim_z920_
#gen100/sim_z1020_
#gen200/sim_z122_g200_
#gen200/sim_z222_g200_
#gen200/sim_z322_g200_
#gen200/sim_z422_g200_
#gen200/sim_z522_g200_
#gen200/sim_z622_g200_
#gen200/sim_z722_g200_
#gen200/sim_z822_g200_
#gen200/sim_z922_g200_
#gen200/sim_z1022_g200_
