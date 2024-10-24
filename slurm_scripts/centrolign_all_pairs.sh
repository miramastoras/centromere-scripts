#!/bin/bash
#
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory
#
# Job name:
#SBATCH --job-name=jeizenga-pairwise-centrolign
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
#SBATCH --array=1-703
#SBATCH --output=array_job_%A_task_%a.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=1:00:00

date
hostname
pwd

#CENTROLIGN=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/centrolign-94eccde
CENTROLIGN=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/centrolign-f3438dd
CHROMDIR=/private/groups/patenlab/jeizenga/centromere/chr12/
FASTADIR=$CHROMDIR/fastas/
WORKDIR=$CHROMDIR/work/
OUTDIR=$CHROMDIR/direct_pairwise_20231127/

mkdir -p $OUTDIR
mkdir -p $WORKDIR
cd $WORKDIR

SAMPLE1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/sample_pair.in_tree.txt | cut -f1)
SAMPLE2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/sample_pair.in_tree.txt | cut -f2)
echo "sample 1:" $SAMPLE1
echo "sample 2:" $SAMPLE2
echo "out:" $OUTDIR

IDENTIFIER1=${SAMPLE1%.*}
IDENTIFIER2=${SAMPLE2%.*}

if [ ${SAMPLE1##*.} == 1 ];
then
    PARENT1=maternal
else
    PARENT1=paternal
fi

if [ ${SAMPLE2##*.} == 1 ];
then
    PARENT2=maternal
else
    PARENT2=paternal
fi

FASTA1=${FASTADIR}/${IDENTIFIER1}.${PARENT1}.hor_array.fasta
FASTA2=${FASTADIR}/${IDENTIFIER2}.${PARENT2}.hor_array.fasta

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2
TEMP_FASTA=${WORKDIR}/${SAMPLE1}_${SAMPLE2}.fa
cat $FASTA1 $FASTA2 > $TEMP_FASTA
${CENTROLIGN} -v 3 --skip-calibration $TEMP_FASTA > $OUTDIR/cigar.${SAMPLE1}_${SAMPLE2}.txt
rm $TEMP_FASTA

#contents of sample_pair.in_tree.txt
#HG00438.1	HG00438.2
#HG00438.1	HG00735.1
#HG00438.1	HG00735.2
#HG00438.1	HG00741.1
#HG00438.1	HG00741.2
#HG00438.1	HG01071.1
#HG00438.1	HG01071.2
#HG00438.1	HG01106.2
#HG00438.1	HG01109.1
#HG00438.1	HG01109.2
#HG00438.1	HG01243.2
#HG00438.1	HG01258.1
#HG00438.1	HG01358.1
#HG00438.1	HG01358.2
#HG00438.1	HG01361.1
#HG00438.1	HG01928.1
#HG00438.1	HG01928.2
#HG00438.1	HG01952.1
#HG00438.1	HG01978.1
#HG00438.1	HG02055.1
#HG00438.1	HG02145.1
#HG00438.1	HG02572.2
#HG00438.1	HG02622.2
#HG00438.1	HG02630.1
#HG00438.1	HG02630.2
#HG00438.1	HG02723.1
#HG00438.1	HG02723.2
#HG00438.1	HG02818.2
#HG00438.1	HG03098.1
#HG00438.1	HG03486.1
#HG00438.1	HG03486.2
#HG00438.1	HG03492.2
#HG00438.1	HG03540.1
#HG00438.1	HG03540.2
#HG00438.1	HG03579.1
#HG00438.1	NA18906.1
#HG00438.1	NA20129.1
#HG00438.2	HG00735.1
#HG00438.2	HG00735.2
#HG00438.2	HG00741.1
#HG00438.2	HG00741.2
#HG00438.2	HG01071.1
#HG00438.2	HG01071.2
#HG00438.2	HG01106.2
#HG00438.2	HG01109.1
#HG00438.2	HG01109.2
#HG00438.2	HG01243.2
#HG00438.2	HG01258.1
#HG00438.2	HG01358.1
#HG00438.2	HG01358.2
#HG00438.2	HG01361.1
#HG00438.2	HG01928.1
#HG00438.2	HG01928.2
#HG00438.2	HG01952.1
#HG00438.2	HG01978.1
#HG00438.2	HG02055.1
#HG00438.2	HG02145.1
#HG00438.2	HG02572.2
#HG00438.2	HG02622.2
#HG00438.2	HG02630.1
#HG00438.2	HG02630.2
#HG00438.2	HG02723.1
#HG00438.2	HG02723.2
#HG00438.2	HG02818.2
#HG00438.2	HG03098.1
#HG00438.2	HG03486.1
#HG00438.2	HG03486.2
#HG00438.2	HG03492.2
#HG00438.2	HG03540.1
#HG00438.2	HG03540.2
#HG00438.2	HG03579.1
#HG00438.2	NA18906.1
#HG00438.2	NA20129.1
#HG00735.1	HG00735.2
#HG00735.1	HG00741.1
#HG00735.1	HG00741.2
#HG00735.1	HG01071.1
#HG00735.1	HG01071.2
#HG00735.1	HG01106.2
#HG00735.1	HG01109.1
#HG00735.1	HG01109.2
#HG00735.1	HG01243.2
#HG00735.1	HG01258.1
#HG00735.1	HG01358.1
#HG00735.1	HG01358.2
#HG00735.1	HG01361.1
#HG00735.1	HG01928.1
#HG00735.1	HG01928.2
#HG00735.1	HG01952.1
#HG00735.1	HG01978.1
#HG00735.1	HG02055.1
#HG00735.1	HG02145.1
#HG00735.1	HG02572.2
#HG00735.1	HG02622.2
#HG00735.1	HG02630.1
#HG00735.1	HG02630.2
#HG00735.1	HG02723.1
#HG00735.1	HG02723.2
#HG00735.1	HG02818.2
#HG00735.1	HG03098.1
#HG00735.1	HG03486.1
#HG00735.1	HG03486.2
#HG00735.1	HG03492.2
#HG00735.1	HG03540.1
#HG00735.1	HG03540.2
#HG00735.1	HG03579.1
#HG00735.1	NA18906.1
#HG00735.1	NA20129.1
#HG00735.2	HG00741.1
#HG00735.2	HG00741.2
#HG00735.2	HG01071.1
#HG00735.2	HG01071.2
#HG00735.2	HG01106.2
#HG00735.2	HG01109.1
#HG00735.2	HG01109.2
#HG00735.2	HG01243.2
#HG00735.2	HG01258.1
#HG00735.2	HG01358.1
#HG00735.2	HG01358.2
#HG00735.2	HG01361.1
#HG00735.2	HG01928.1
#HG00735.2	HG01928.2
#HG00735.2	HG01952.1
#HG00735.2	HG01978.1
#HG00735.2	HG02055.1
#HG00735.2	HG02145.1
#HG00735.2	HG02572.2
#HG00735.2	HG02622.2
#HG00735.2	HG02630.1
#HG00735.2	HG02630.2
#HG00735.2	HG02723.1
#HG00735.2	HG02723.2
#HG00735.2	HG02818.2
#HG00735.2	HG03098.1
#HG00735.2	HG03486.1
#HG00735.2	HG03486.2
#HG00735.2	HG03492.2
#HG00735.2	HG03540.1
#HG00735.2	HG03540.2
#HG00735.2	HG03579.1
#HG00735.2	NA18906.1
#HG00735.2	NA20129.1
#HG00741.1	HG00741.2
#HG00741.1	HG01071.1
#HG00741.1	HG01071.2
#HG00741.1	HG01106.2
#HG00741.1	HG01109.1
#HG00741.1	HG01109.2
#HG00741.1	HG01243.2
#HG00741.1	HG01258.1
#HG00741.1	HG01358.1
#HG00741.1	HG01358.2
#HG00741.1	HG01361.1
#HG00741.1	HG01928.1
#HG00741.1	HG01928.2
#HG00741.1	HG01952.1
#HG00741.1	HG01978.1
#HG00741.1	HG02055.1
#HG00741.1	HG02145.1
#HG00741.1	HG02572.2
#HG00741.1	HG02622.2
#HG00741.1	HG02630.1
#HG00741.1	HG02630.2
#HG00741.1	HG02723.1
#HG00741.1	HG02723.2
#HG00741.1	HG02818.2
#HG00741.1	HG03098.1
#HG00741.1	HG03486.1
#HG00741.1	HG03486.2
#HG00741.1	HG03492.2
#HG00741.1	HG03540.1
#HG00741.1	HG03540.2
#HG00741.1	HG03579.1
#HG00741.1	NA18906.1
#HG00741.1	NA20129.1
#HG00741.2	HG01071.1
#HG00741.2	HG01071.2
#HG00741.2	HG01106.2
#HG00741.2	HG01109.1
#HG00741.2	HG01109.2
#HG00741.2	HG01243.2
#HG00741.2	HG01258.1
#HG00741.2	HG01358.1
#HG00741.2	HG01358.2
#HG00741.2	HG01361.1
#HG00741.2	HG01928.1
#HG00741.2	HG01928.2
#HG00741.2	HG01952.1
#HG00741.2	HG01978.1
#HG00741.2	HG02055.1
#HG00741.2	HG02145.1
#HG00741.2	HG02572.2
#HG00741.2	HG02622.2
#HG00741.2	HG02630.1
#HG00741.2	HG02630.2
#HG00741.2	HG02723.1
#HG00741.2	HG02723.2
#HG00741.2	HG02818.2
#HG00741.2	HG03098.1
#HG00741.2	HG03486.1
#HG00741.2	HG03486.2
#HG00741.2	HG03492.2
#HG00741.2	HG03540.1
#HG00741.2	HG03540.2
#HG00741.2	HG03579.1
#HG00741.2	NA18906.1
#HG00741.2	NA20129.1
#HG01071.1	HG01071.2
#HG01071.1	HG01106.2
#HG01071.1	HG01109.1
#HG01071.1	HG01109.2
#HG01071.1	HG01243.2
#HG01071.1	HG01258.1
#HG01071.1	HG01358.1
#HG01071.1	HG01358.2
#HG01071.1	HG01361.1
#HG01071.1	HG01928.1
#HG01071.1	HG01928.2
#HG01071.1	HG01952.1
#HG01071.1	HG01978.1
#HG01071.1	HG02055.1
#HG01071.1	HG02145.1
#HG01071.1	HG02572.2
#HG01071.1	HG02622.2
#HG01071.1	HG02630.1
#HG01071.1	HG02630.2
#HG01071.1	HG02723.1
#HG01071.1	HG02723.2
#HG01071.1	HG02818.2
#HG01071.1	HG03098.1
#HG01071.1	HG03486.1
#HG01071.1	HG03486.2
#HG01071.1	HG03492.2
#HG01071.1	HG03540.1
#HG01071.1	HG03540.2
#HG01071.1	HG03579.1
#HG01071.1	NA18906.1
#HG01071.1	NA20129.1
#HG01071.2	HG01106.2
#HG01071.2	HG01109.1
#HG01071.2	HG01109.2
#HG01071.2	HG01243.2
#HG01071.2	HG01258.1
#HG01071.2	HG01358.1
#HG01071.2	HG01358.2
#HG01071.2	HG01361.1
#HG01071.2	HG01928.1
#HG01071.2	HG01928.2
#HG01071.2	HG01952.1
#HG01071.2	HG01978.1
#HG01071.2	HG02055.1
#HG01071.2	HG02145.1
#HG01071.2	HG02572.2
#HG01071.2	HG02622.2
#HG01071.2	HG02630.1
#HG01071.2	HG02630.2
#HG01071.2	HG02723.1
#HG01071.2	HG02723.2
#HG01071.2	HG02818.2
#HG01071.2	HG03098.1
#HG01071.2	HG03486.1
#HG01071.2	HG03486.2
#HG01071.2	HG03492.2
#HG01071.2	HG03540.1
#HG01071.2	HG03540.2
#HG01071.2	HG03579.1
#HG01071.2	NA18906.1
#HG01071.2	NA20129.1
#HG01106.2	HG01109.1
#HG01106.2	HG01109.2
#HG01106.2	HG01243.2
#HG01106.2	HG01258.1
#HG01106.2	HG01358.1
#HG01106.2	HG01358.2
#HG01106.2	HG01361.1
#HG01106.2	HG01928.1
#HG01106.2	HG01928.2
#HG01106.2	HG01952.1
#HG01106.2	HG01978.1
#HG01106.2	HG02055.1
#HG01106.2	HG02145.1
#HG01106.2	HG02572.2
#HG01106.2	HG02622.2
#HG01106.2	HG02630.1
#HG01106.2	HG02630.2
#HG01106.2	HG02723.1
#HG01106.2	HG02723.2
#HG01106.2	HG02818.2
#HG01106.2	HG03098.1
#HG01106.2	HG03486.1
#HG01106.2	HG03486.2
#HG01106.2	HG03492.2
#HG01106.2	HG03540.1
#HG01106.2	HG03540.2
#HG01106.2	HG03579.1
#HG01106.2	NA18906.1
#HG01106.2	NA20129.1
#HG01109.1	HG01109.2
#HG01109.1	HG01243.2
#HG01109.1	HG01258.1
#HG01109.1	HG01358.1
#HG01109.1	HG01358.2
#HG01109.1	HG01361.1
#HG01109.1	HG01928.1
#HG01109.1	HG01928.2
#HG01109.1	HG01952.1
#HG01109.1	HG01978.1
#HG01109.1	HG02055.1
#HG01109.1	HG02145.1
#HG01109.1	HG02572.2
#HG01109.1	HG02622.2
#HG01109.1	HG02630.1
#HG01109.1	HG02630.2
#HG01109.1	HG02723.1
#HG01109.1	HG02723.2
#HG01109.1	HG02818.2
#HG01109.1	HG03098.1
#HG01109.1	HG03486.1
#HG01109.1	HG03486.2
#HG01109.1	HG03492.2
#HG01109.1	HG03540.1
#HG01109.1	HG03540.2
#HG01109.1	HG03579.1
#HG01109.1	NA18906.1
#HG01109.1	NA20129.1
#HG01109.2	HG01243.2
#HG01109.2	HG01258.1
#HG01109.2	HG01358.1
#HG01109.2	HG01358.2
#HG01109.2	HG01361.1
#HG01109.2	HG01928.1
#HG01109.2	HG01928.2
#HG01109.2	HG01952.1
#HG01109.2	HG01978.1
#HG01109.2	HG02055.1
#HG01109.2	HG02145.1
#HG01109.2	HG02572.2
#HG01109.2	HG02622.2
#HG01109.2	HG02630.1
#HG01109.2	HG02630.2
#HG01109.2	HG02723.1
#HG01109.2	HG02723.2
#HG01109.2	HG02818.2
#HG01109.2	HG03098.1
#HG01109.2	HG03486.1
#HG01109.2	HG03486.2
#HG01109.2	HG03492.2
#HG01109.2	HG03540.1
#HG01109.2	HG03540.2
#HG01109.2	HG03579.1
#HG01109.2	NA18906.1
#HG01109.2	NA20129.1
#HG01243.2	HG01258.1
#HG01243.2	HG01358.1
#HG01243.2	HG01358.2
#HG01243.2	HG01361.1
#HG01243.2	HG01928.1
#HG01243.2	HG01928.2
#HG01243.2	HG01952.1
#HG01243.2	HG01978.1
#HG01243.2	HG02055.1
#HG01243.2	HG02145.1
#HG01243.2	HG02572.2
#HG01243.2	HG02622.2
#HG01243.2	HG02630.1
#HG01243.2	HG02630.2
#HG01243.2	HG02723.1
#HG01243.2	HG02723.2
#HG01243.2	HG02818.2
#HG01243.2	HG03098.1
#HG01243.2	HG03486.1
#HG01243.2	HG03486.2
#HG01243.2	HG03492.2
#HG01243.2	HG03540.1
#HG01243.2	HG03540.2
#HG01243.2	HG03579.1
#HG01243.2	NA18906.1
#HG01243.2	NA20129.1
#HG01258.1	HG01358.1
#HG01258.1	HG01358.2
#HG01258.1	HG01361.1
#HG01258.1	HG01928.1
#HG01258.1	HG01928.2
#HG01258.1	HG01952.1
#HG01258.1	HG01978.1
#HG01258.1	HG02055.1
#HG01258.1	HG02145.1
#HG01258.1	HG02572.2
#HG01258.1	HG02622.2
#HG01258.1	HG02630.1
#HG01258.1	HG02630.2
#HG01258.1	HG02723.1
#HG01258.1	HG02723.2
#HG01258.1	HG02818.2
#HG01258.1	HG03098.1
#HG01258.1	HG03486.1
#HG01258.1	HG03486.2
#HG01258.1	HG03492.2
#HG01258.1	HG03540.1
#HG01258.1	HG03540.2
#HG01258.1	HG03579.1
#HG01258.1	NA18906.1
#HG01258.1	NA20129.1
#HG01358.1	HG01358.2
#HG01358.1	HG01361.1
#HG01358.1	HG01928.1
#HG01358.1	HG01928.2
#HG01358.1	HG01952.1
#HG01358.1	HG01978.1
#HG01358.1	HG02055.1
#HG01358.1	HG02145.1
#HG01358.1	HG02572.2
#HG01358.1	HG02622.2
#HG01358.1	HG02630.1
#HG01358.1	HG02630.2
#HG01358.1	HG02723.1
#HG01358.1	HG02723.2
#HG01358.1	HG02818.2
#HG01358.1	HG03098.1
#HG01358.1	HG03486.1
#HG01358.1	HG03486.2
#HG01358.1	HG03492.2
#HG01358.1	HG03540.1
#HG01358.1	HG03540.2
#HG01358.1	HG03579.1
#HG01358.1	NA18906.1
#HG01358.1	NA20129.1
#HG01358.2	HG01361.1
#HG01358.2	HG01928.1
#HG01358.2	HG01928.2
#HG01358.2	HG01952.1
#HG01358.2	HG01978.1
#HG01358.2	HG02055.1
#HG01358.2	HG02145.1
#HG01358.2	HG02572.2
#HG01358.2	HG02622.2
#HG01358.2	HG02630.1
#HG01358.2	HG02630.2
#HG01358.2	HG02723.1
#HG01358.2	HG02723.2
#HG01358.2	HG02818.2
#HG01358.2	HG03098.1
#HG01358.2	HG03486.1
#HG01358.2	HG03486.2
#HG01358.2	HG03492.2
#HG01358.2	HG03540.1
#HG01358.2	HG03540.2
#HG01358.2	HG03579.1
#HG01358.2	NA18906.1
#HG01358.2	NA20129.1
#HG01361.1	HG01928.1
#HG01361.1	HG01928.2
#HG01361.1	HG01952.1
#HG01361.1	HG01978.1
#HG01361.1	HG02055.1
#HG01361.1	HG02145.1
#HG01361.1	HG02572.2
#HG01361.1	HG02622.2
#HG01361.1	HG02630.1
#HG01361.1	HG02630.2
#HG01361.1	HG02723.1
#HG01361.1	HG02723.2
#HG01361.1	HG02818.2
#HG01361.1	HG03098.1
#HG01361.1	HG03486.1
#HG01361.1	HG03486.2
#HG01361.1	HG03492.2
#HG01361.1	HG03540.1
#HG01361.1	HG03540.2
#HG01361.1	HG03579.1
#HG01361.1	NA18906.1
#HG01361.1	NA20129.1
#HG01928.1	HG01928.2
#HG01928.1	HG01952.1
#HG01928.1	HG01978.1
#HG01928.1	HG02055.1
#HG01928.1	HG02145.1
#HG01928.1	HG02572.2
#HG01928.1	HG02622.2
#HG01928.1	HG02630.1
#HG01928.1	HG02630.2
#HG01928.1	HG02723.1
#HG01928.1	HG02723.2
#HG01928.1	HG02818.2
#HG01928.1	HG03098.1
#HG01928.1	HG03486.1
#HG01928.1	HG03486.2
#HG01928.1	HG03492.2
#HG01928.1	HG03540.1
#HG01928.1	HG03540.2
#HG01928.1	HG03579.1
#HG01928.1	NA18906.1
#HG01928.1	NA20129.1
#HG01928.2	HG01952.1
#HG01928.2	HG01978.1
#HG01928.2	HG02055.1
#HG01928.2	HG02145.1
#HG01928.2	HG02572.2
#HG01928.2	HG02622.2
#HG01928.2	HG02630.1
#HG01928.2	HG02630.2
#HG01928.2	HG02723.1
#HG01928.2	HG02723.2
#HG01928.2	HG02818.2
#HG01928.2	HG03098.1
#HG01928.2	HG03486.1
#HG01928.2	HG03486.2
#HG01928.2	HG03492.2
#HG01928.2	HG03540.1
#HG01928.2	HG03540.2
#HG01928.2	HG03579.1
#HG01928.2	NA18906.1
#HG01928.2	NA20129.1
#HG01952.1	HG01978.1
#HG01952.1	HG02055.1
#HG01952.1	HG02145.1
#HG01952.1	HG02572.2
#HG01952.1	HG02622.2
#HG01952.1	HG02630.1
#HG01952.1	HG02630.2
#HG01952.1	HG02723.1
#HG01952.1	HG02723.2
#HG01952.1	HG02818.2
#HG01952.1	HG03098.1
#HG01952.1	HG03486.1
#HG01952.1	HG03486.2
#HG01952.1	HG03492.2
#HG01952.1	HG03540.1
#HG01952.1	HG03540.2
#HG01952.1	HG03579.1
#HG01952.1	NA18906.1
#HG01952.1	NA20129.1
#HG01978.1	HG02055.1
#HG01978.1	HG02145.1
#HG01978.1	HG02572.2
#HG01978.1	HG02622.2
#HG01978.1	HG02630.1
#HG01978.1	HG02630.2
#HG01978.1	HG02723.1
#HG01978.1	HG02723.2
#HG01978.1	HG02818.2
#HG01978.1	HG03098.1
#HG01978.1	HG03486.1
#HG01978.1	HG03486.2
#HG01978.1	HG03492.2
#HG01978.1	HG03540.1
#HG01978.1	HG03540.2
#HG01978.1	HG03579.1
#HG01978.1	NA18906.1
#HG01978.1	NA20129.1
#HG02055.1	HG02145.1
#HG02055.1	HG02572.2
#HG02055.1	HG02622.2
#HG02055.1	HG02630.1
#HG02055.1	HG02630.2
#HG02055.1	HG02723.1
#HG02055.1	HG02723.2
#HG02055.1	HG02818.2
#HG02055.1	HG03098.1
#HG02055.1	HG03486.1
#HG02055.1	HG03486.2
#HG02055.1	HG03492.2
#HG02055.1	HG03540.1
#HG02055.1	HG03540.2
#HG02055.1	HG03579.1
#HG02055.1	NA18906.1
#HG02055.1	NA20129.1
#HG02145.1	HG02572.2
#HG02145.1	HG02622.2
#HG02145.1	HG02630.1
#HG02145.1	HG02630.2
#HG02145.1	HG02723.1
#HG02145.1	HG02723.2
#HG02145.1	HG02818.2
#HG02145.1	HG03098.1
#HG02145.1	HG03486.1
#HG02145.1	HG03486.2
#HG02145.1	HG03492.2
#HG02145.1	HG03540.1
#HG02145.1	HG03540.2
#HG02145.1	HG03579.1
#HG02145.1	NA18906.1
#HG02145.1	NA20129.1
#HG02572.2	HG02622.2
#HG02572.2	HG02630.1
#HG02572.2	HG02630.2
#HG02572.2	HG02723.1
#HG02572.2	HG02723.2
#HG02572.2	HG02818.2
#HG02572.2	HG03098.1
#HG02572.2	HG03486.1
#HG02572.2	HG03486.2
#HG02572.2	HG03492.2
#HG02572.2	HG03540.1
#HG02572.2	HG03540.2
#HG02572.2	HG03579.1
#HG02572.2	NA18906.1
#HG02572.2	NA20129.1
#HG02622.2	HG02630.1
#HG02622.2	HG02630.2
#HG02622.2	HG02723.1
#HG02622.2	HG02723.2
#HG02622.2	HG02818.2
#HG02622.2	HG03098.1
#HG02622.2	HG03486.1
#HG02622.2	HG03486.2
#HG02622.2	HG03492.2
#HG02622.2	HG03540.1
#HG02622.2	HG03540.2
#HG02622.2	HG03579.1
#HG02622.2	NA18906.1
#HG02622.2	NA20129.1
#HG02630.1	HG02630.2
#HG02630.1	HG02723.1
#HG02630.1	HG02723.2
#HG02630.1	HG02818.2
#HG02630.1	HG03098.1
#HG02630.1	HG03486.1
#HG02630.1	HG03486.2
#HG02630.1	HG03492.2
#HG02630.1	HG03540.1
#HG02630.1	HG03540.2
#HG02630.1	HG03579.1
#HG02630.1	NA18906.1
#HG02630.1	NA20129.1
#HG02630.2	HG02723.1
#HG02630.2	HG02723.2
#HG02630.2	HG02818.2
#HG02630.2	HG03098.1
#HG02630.2	HG03486.1
#HG02630.2	HG03486.2
#HG02630.2	HG03492.2
#HG02630.2	HG03540.1
#HG02630.2	HG03540.2
#HG02630.2	HG03579.1
#HG02630.2	NA18906.1
#HG02630.2	NA20129.1
#HG02723.1	HG02723.2
#HG02723.1	HG02818.2
#HG02723.1	HG03098.1
#HG02723.1	HG03486.1
#HG02723.1	HG03486.2
#HG02723.1	HG03492.2
#HG02723.1	HG03540.1
#HG02723.1	HG03540.2
#HG02723.1	HG03579.1
#HG02723.1	NA18906.1
#HG02723.1	NA20129.1
#HG02723.2	HG02818.2
#HG02723.2	HG03098.1
#HG02723.2	HG03486.1
#HG02723.2	HG03486.2
#HG02723.2	HG03492.2
#HG02723.2	HG03540.1
#HG02723.2	HG03540.2
#HG02723.2	HG03579.1
#HG02723.2	NA18906.1
#HG02723.2	NA20129.1
#HG02818.2	HG03098.1
#HG02818.2	HG03486.1
#HG02818.2	HG03486.2
#HG02818.2	HG03492.2
#HG02818.2	HG03540.1
#HG02818.2	HG03540.2
#HG02818.2	HG03579.1
#HG02818.2	NA18906.1
#HG02818.2	NA20129.1
#HG03098.1	HG03486.1
#HG03098.1	HG03486.2
#HG03098.1	HG03492.2
#HG03098.1	HG03540.1
#HG03098.1	HG03540.2
#HG03098.1	HG03579.1
#HG03098.1	NA18906.1
#HG03098.1	NA20129.1
#HG03486.1	HG03486.2
#HG03486.1	HG03492.2
#HG03486.1	HG03540.1
#HG03486.1	HG03540.2
#HG03486.1	HG03579.1
#HG03486.1	NA18906.1
#HG03486.1	NA20129.1
#HG03486.2	HG03492.2
#HG03486.2	HG03540.1
#HG03486.2	HG03540.2
#HG03486.2	HG03579.1
#HG03486.2	NA18906.1
#HG03486.2	NA20129.1
#HG03492.2	HG03540.1
#HG03492.2	HG03540.2
#HG03492.2	HG03579.1
#HG03492.2	NA18906.1
#HG03492.2	NA20129.1
#HG03540.1	HG03540.2
#HG03540.1	HG03579.1
#HG03540.1	NA18906.1
#HG03540.1	NA20129.1
#HG03540.2	HG03579.1
#HG03540.2	NA18906.1
#HG03540.2	NA20129.1
#HG03579.1	NA18906.1
#HG03579.1	NA20129.1
#NA18906.1	NA20129.1