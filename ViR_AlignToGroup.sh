#!/bin/bash -x

#----------------------------------------------------------------------------------------------------------------------------
# Input settings
#----------------------------------------------------------------------------------------------------------------------------

while [[ $# > 0 ]]
do
key="$1"

case $key in
    # sample sam = Reads_Alignments.txt
    -sample_sam)
    SAMPLE_SAM="$2"
    shift # past argument
    ;;
    -sample_name)
    SAMPLE_NAME="$2"
    shift # past argument
    ;;
    -work_files_dir)
    WORKDIR="$2"
    shift # past argument
    ;;
    -group_name)
    GROUP="$2"
    shift
    ;;
    -groups_fasta)
    GROUPS_FASTA="$2"
    shift
    ;;
    -reads_list)
    READ_LIST="$2"
    shift # past argument
    ;;
    -outdir)
    OUTDIR="$2"
    shift # past argument
    ;;
    -bedtools_exe)
    BEDTOOLS="$2"
    shift # past argument
    ;;
    -trinity_exe)
    TRINITY="$2"
    shift # past argument
    ;;
    -samtools_exe)
    SAMTOOLS="$2"
    shift # past argument
    ;;
    -bwa_exe)
    BWA="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
    # unknown option
    ;;

esac
shift # past argument or value
done

#----------------------------------------------------------------------------------------------------------------------------
# Concatenate sam file in case of several lanes for the same sample
#----------------------------------------------------------------------------------------------------------------------------
python $WORKDIR/RC_HeaderSam.py -sam $SAMPLE_SAM > $OUTDIR/${GROUP}_${SAMPLE_NAME}.sam

while IFS='' read -r line || [[ -n "$line" ]]; do
	echo -e $line
	grep -w $line $SAMPLE_SAM >> $OUTDIR/${GROUP}_${SAMPLE_NAME}.sam
	
done < $READ_LIST

#----------------------------------------------------------------------------------------------------------------------------
# Converting Sam to Bam
#----------------------------------------------------------------------------------------------------------------------------

$SAMTOOLS view -bS $OUTDIR/${GROUP}_${SAMPLE_NAME}.sam > $OUTDIR/${GROUP}_${SAMPLE_NAME}.bam

#----------------------------------------------------------------------------------------------------------------------------
# Sorting Bam
#----------------------------------------------------------------------------------------------------------------------------

$SAMTOOLS sort -n -o $OUTDIR/${GROUP}_${SAMPLE_NAME}_sorted.bam $OUTDIR/${GROUP}_${SAMPLE_NAME}.bam

#----------------------------------------------------------------------------------------------------------------------------
# Fastq editing
#----------------------------------------------------------------------------------------------------------------------------

$BEDTOOLS bamtofastq \
-i $OUTDIR/${GROUP}_${SAMPLE_NAME}_sorted.bam \
-fq $OUTDIR/${GROUP}_${SAMPLE_NAME}_end1.fq \
-fq2 $OUTDIR/${GROUP}_${SAMPLE_NAME}_end2.fq > $OUTDIR/${GROUP}_${SAMPLE_NAME}_stdout.txt

#----------------------------------------------------------------------------------------------------------------------------
# Alingment Reads to the equivalent region
#----------------------------------------------------------------------------------------------------------------------------
$BWA index $GROUPS_FASTA >> $OUTDIR/${GROUP}_${SAMPLE_NAME}_stdout.txt
$BWA mem $GROUPS_FASTA $OUTDIR/${GROUP}_${SAMPLE_NAME}_end1.fq $OUTDIR/${GROUP}_${SAMPLE_NAME}_end2.fq > $OUTDIR/${GROUP}_${SAMPLE_NAME}_secAlign.sam
$SAMTOOLS view -Sb  $OUTDIR/${GROUP}_${SAMPLE_NAME}_secAlign.sam >  $OUTDIR/${GROUP}_${SAMPLE_NAME}_secAlign.bam
$SAMTOOLS sort -o $OUTDIR/${GROUP}_${SAMPLE_NAME}_secAlign_sortedIGV.bam $OUTDIR/${GROUP}_${SAMPLE_NAME}_secAlign.bam >> $OUTDIR/${GROUP}_${SAMPLE_NAME}_stdout.txt
$SAMTOOLS index $OUTDIR/${GROUP}_${SAMPLE_NAME}_secAlign_sortedIGV.bam >> $OUTDIR/${GROUP}_${SAMPLE_NAME}_stdout.txt

#----------------------------------------------------------------------------------------------------------------------------
# Try de novo Assembly
#----------------------------------------------------------------------------------------------------------------------------
$TRINITY --seqType fq \
--max_memory 10G \
--left $OUTDIR/${GROUP}_${SAMPLE_NAME}_end1.fq \
--right $OUTDIR/${GROUP}_${SAMPLE_NAME}_end2.fq \
--output $OUTDIR/trinity \
--full_cleanup >> $OUTDIR/${GROUP}_${SAMPLE_NAME}_stdout.txt