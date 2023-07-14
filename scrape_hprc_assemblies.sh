#!/bin/bash
# Job name:
#SBATCH --job-name=jeizenga-scrape-active-arrays
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
#SBATCH --array=1-47
#SBATCH --output=array_job_%A_task_%a.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=4:00:00
#
## script to run:

date
hostname
pwd

EXTRACTDIR=/private/groups/patenlab/jeizenga/centromere/extraction/

# get the sample ID from the config
SAMPLE=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$EXTRACTDIR"/samples.txt)
echo "sample:" $SAMPLE

# the web location of the assemblies
ASSEMBLY_ENDPOINT="https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/728E4476-8D84-4E8E-BA6D-AC9BF482ECCD--YEAR_1_GENBANK_ASSEM/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank"

# the wb location of the HOR annotations
ANNOTATION_ENDPOINT="https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/08934468-0AE3-42B6-814A-C5422311A53D--HUMAS_HMMER/$SAMPLE"

echo "assembly endpoint:" $ASSEMBLY_ENDPOINT
echo "annotation endpoint:" $ANNOTATION_ENDPOINT

cd "$EXTRACTDIR"/work/

for PARENT in maternal paternal;
do
    ASSEMBLY="$SAMPLE"."$PARENT".f1_assembly_v2_genbank.fa
    HUMAS_BED=AS-HOR-vs-"$SAMPLE"-"$PARENT".bed
    
    OUTDIR="$EXTRACTDIR"/outputs/beds/"$SAMPLE"/
    mkdir -p "$OUTDIR"
    
    # download and extract data
    echo "downloading assembly" $ASSEMBLY
    wget -q "$ASSEMBLY_ENDPOINT"/"$ASSEMBLY".gz
    echo "downloading annotations" $HUMAS_BED
    wget -q "$ANNOTATION_ENDPOINT"/"$HUMAS_BED"
    echo "decompressing assembly"
    gunzip "$ASSEMBLY".gz

    # find the HORs
    echo "looking for HORs for" $PARENT "haplotype of" $SAMPLE
    /private/groups/patenlab/jeizenga/GitHub/centromere-scripts/locate_hor_arrays.py \
        -a "$ASSEMBLY" \
        -f "$EXTRACTDIR"/new_flanks_unique_200000.fasta \
        -b "$HUMAS_BED" \
        > "$OUTDIR"/hor_arrays."$PARENT".bed
    
    # clean up the files
    echo "cleaning up data"
    rm "$ASSEMBLY" "$HUMAS_BED"
done

# contents of samples.txt:
#HG002
#HG00438
#HG005
#HG00621
#HG00673
#HG00733
#HG00735
#HG00741
#HG01071
#HG01106
#HG01109
#HG01123
#HG01175
#HG01243
#HG01258
#HG01358
#HG01361
#HG01891
#HG01928
#HG01952
#HG01978
#HG02055
#HG02080
#HG02109
#HG02145
#HG02148
#HG02257
#HG02486
#HG02559
#HG02572
#HG02622
#HG02630
#HG02717
#HG02723
#HG02818
#HG02886
#HG03098
#HG03453
#HG03486
#HG03492
#HG03516
#HG03540
#HG03579
#NA18906
#NA19240
#NA20129
#NA21309
