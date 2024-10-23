#!/bin/bash
#
# Slurm script intended to produce pairwise alignments of the WashU pedigree centromeres using centrolign
#
# Job name:
#SBATCH --job-name=jeizenga-pedigree-centrolign
#
# Partition - This is the queue it goes in:
#SBATCH --partition=short
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
#SBATCH --array=1-211
#SBATCH --output=array_job_%A_task_%a.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=1:00:00

date
hostname
pwd


PEDDIR=/private/groups/patenlab/jeizenga/centromere/washu_noduplex/
PAIRS=$PEDDIR/paired_alns.txt

OUTDIR=$PEDDIR/alignments/
mkdir -p $OUTDIR

CHILD_ARRAY=$(awk "NR==$SLURM_ARRAY_TASK_ID" $PAIRS | cut -f 1)
RELATIVE_ARRAY=$(awk "NR==$SLURM_ARRAY_TASK_ID" $PAIRS | cut -f 2)
CHILD_FASTA=$PEDDIR/$CHILD_ARRAY
RELATIVE_FASTA=$PEDDIR/$RELATIVE_ARRAY
RELATIVE=$(echo $RELATIVE_ARRAY | grep -Eo "PAN[0-9]+")
RELATIVE_HAP=$(echo $RELATIVE_ARRAY | grep -Eo "haplotype[0-9]")
CHILD=$(echo $CHILD_ARRAY | grep -Eo "PAN[0-9]+")
CHILD_HAP=$(echo $CHILD_ARRAY | grep -Eo "haplotype[0-9]")
CHR=$(echo $CHILD_ARRAY | grep -Eo "chr[^\.]+")

OUTFILE=$OUTDIR/"$CHILD"_"$CHILD_HAP"_"$RELATIVE"_"$RELATIVE_HAP"."$CHR".cigar.txt

echo "input"
echo $CHILD_FASTA
echo $RELATIVE_FASTA
echo "relative" $RELATIVE "haplotype" $RELATIVE_HAP
echo "child" $CHILD "haplotype" $CHILD_HAP
echo "output:"
echo $OUTFILE

CENTROLIGN_DIR=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/
CENTROLIGN=$CENTROLIGN_DIR/centrolign-86d3252

WORKDIR=$PEDDIR/work/
mkdir -p $WORKDIR
cd $WORKDIR

echo "aligning with centrolign"
TEMP_FASTA=${WORKDIR}/joined_"$SLURM_ARRAY_TASK_ID".fa
mkdir -p `dirname $TEMP_FASTA`
cat $CHILD_FASTA $RELATIVE_FASTA > $TEMP_FASTA
mkdir -p `dirname $OUTFILE`
/usr/bin/time -v ${CENTROLIGN} -v 3 $TEMP_FASTA > $OUTFILE
rm $TEMP_FASTA


