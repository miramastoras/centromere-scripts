#!/bin/bash
#
# Slurm script to download HPRC assemblies and extract the centromeres
# as provided by a BED file
#
# Job name:
#SBATCH --job-name=jeizenga-extract-hor-arrays
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
#SBATCH --mem=1gb
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

CHR=chr12

EXTRACTDIR=/private/groups/patenlab/jeizenga/centromere/extraction/
WORKDIR=$EXTRACTDIR/work/
BEDDIR=$EXTRACTDIR/outputs/beds/
OUTDIR=$EXTRACTDIR/outputs/hor_arrays/$CHR/

mkdir -p $OUTDIR
cd $WORKDIR

SAMPLE=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$EXTRACTDIR"/samples.txt)
echo "sample:" $SAMPLE
echo "chrom:" $CHR
echo "out:" $OUTDIR

# the web location of the assemblies
ASSEMBLY_ENDPOINT="https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/728E4476-8D84-4E8E-BA6D-AC9BF482ECCD--YEAR_1_GENBANK_ASSEM/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank"

for PARENT in maternal paternal;
do
    
    ASSEMBLY="$SAMPLE"."$PARENT".f1_assembly_v2_genbank.fa
    BED=$BEDDIR/$SAMPLE/hor_arrays."$PARENT".bed
    
    if [ $PARENT = maternal ];
    then
        PARNUM=1
    else
        PARNUM=2
    fi
    
    cd $WORKDIR
    
    REGIONFILE="$SAMPLE"."$PARENT".txt
    
    grep $CHR $BED | awk '{ printf "%s:%d-%d\n", $1, $2+1, $3 }' > "$REGIONFILE"
    
    STRAND=$(grep $CHR $BED | cut -f 6)
    
    if [ -s $REGIONFILE ];
    then
        
        echo "downloading assembly" $ASSEMBLY
        wget -q "$ASSEMBLY_ENDPOINT"/"$ASSEMBLY".gz
        echo "decompressing assembly"
        gunzip "$ASSEMBLY".gz
                
        # extract and add the sample name as the sequence name
        echo "extract region" `cat $REGIONFILE`
        echo "strand" $STRAND
        
        HORFASTA="$OUTDIR"/"$SAMPLE"."$PARENT".hor_array.fasta
        
        samtools faidx -r $REGIONFILE "$WORKDIR"/"$ASSEMBLY" | sed "s/>/>$SAMPLE.$PARNUM /g" > $HORFASTA
        
        if [ $STRAND = "-" ];
        then
            echo "reverse complementing sequence"
            TMPFILE=$(mktemp)
            /private/groups/patenlab/jeizenga/GitHub/centromere-scripts/fasta_to_rev_comp.py $HORFASTA > $TMPFILE
            mv $TMPFILE $HORFASTA
        fi
        
        rm "$WORKDIR"/"$ASSEMBLY"
        rm "$WORKDIR"/"$ASSEMBLY".fai
    fi
    
    
    rm "$WORKDIR"/"$REGIONFILE"
    
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



