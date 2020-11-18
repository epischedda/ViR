#!/bin/bash -x
#----------------------------------------------------------------------------------------------------------------------------
# Remember to have in path the directories of jellyfish, bowtie e salmon
#----------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------
# Input settings
#----------------------------------------------------------------------------------------------------------------------------

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -th)
    TH="$2"
    shift # past argument
    ;;
    -mem)
    MEM="$2"
    shift # past argument
    ;;
    -working_dir)
    WORKDIR="$2"
    shift # past argument
    ;;
    -non_host_fasta)
    VIRUS="$2"
    shift # past argument
    ;;
    -analysis_name)
    SAMPLE="$2"
    shift # past argument
    ;;
    -read_list)
    SAMPLE_LIST="$2"
    shift # past argument
    ;;
    -trinity_exe)
    TRINITY_EXE="$2"
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

START=$(date +%s);

#----------------------------------------------------------------------------------------------------------------------------
# Control settings
#----------------------------------------------------------------------------------------------------------------------------
if [ -z "$TH" ];then
  TH=1
fi
if [ -z "$MEM" ];then
  MEM=2G
fi

virus_name=$(echo $VIRUS | rev | cut -f 1 -d'/' | rev | cut -f 1 -d'.');
mkdir $WORKDIR/Alignments
mkdir $WORKDIR/BamMerged
mkdir $WORKDIR/Fastq

#----------------------------------------------------------------------------------------------------------------------------
# Index reference sequence (virus.fa)
#----------------------------------------------------------------------------------------------------------------------------
bwa index $VIRUS

#----------------------------------------------------------------------------------------------------------------------------
# For each couple of reads assigned to a sample: 
# Alignment of the sample reads to the virus and extraction of the chimeric reads and couples of reads both mapping to the virus
#----------------------------------------------------------------------------------------------------------------------------

while IFS='' read -r line || [[ -n "$line" ]]; do
    	
	R1=$(echo $line | cut -f1 -d$' ');
	R2=$(echo $line | cut -f2 -d$' ');
	ID_SAMPLE=$(echo $R1 | rev | cut -f 1 -d'/' | rev | cut -f 1 -d'.'| cut -f 1 -d'_');
	LANE=$(echo $R1 | rev | cut -f 1 -d'/' | rev | cut -f 1 -d'.'| cut -f 3 -d'_');
	
	bwa mem -t $TH $VIRUS $R1 $R2 > $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}.sam
	samtools view -@ $TH -Sb  $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}.sam >  $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}.bam
	samtools view -@ $TH -b -F 4 $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}.bam -o $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}_Mapped.bam
	samtools view -@ $TH -u -f 4 -F 264 $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}.bam -o $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}_oneEndUnmapped.bam

	samtools merge -@ $TH -f $WORKDIR/BamMerged/${ID_SAMPLE}_${LANE}_to_${virus_name}_merged.bam $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}_Mapped.bam $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}_oneEndUnmapped.bam
	
	#----------------------------------------------------------------------------------------------------------------------------
	# Remove unused files
	#----------------------------------------------------------------------------------------------------------------------------
	rm $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}.sam
	rm $WORKDIR/Alignments/${ID_SAMPLE}_${LANE}_to_${virus_name}.bam
	
	#----------------------------------------------------------------------------------------------------------------------------
	# Sorting by name of the merged file
	#----------------------------------------------------------------------------------------------------------------------------
	samtools sort -@ $TH -n -o $WORKDIR/BamMerged/${ID_SAMPLE}_${LANE}_to_${virus_name}_sortedName.bam $WORKDIR/BamMerged/${ID_SAMPLE}_${LANE}_to_${virus_name}_merged.bam
	
	#----------------------------------------------------------------------------------------------------------------------------
	# Remove unused files
	#----------------------------------------------------------------------------------------------------------------------------
	rm $WORKDIR/BamMerged/${ID_SAMPLE}_${LANE}_to_${virus_name}_merged.bam

done < $SAMPLE_LIST

#----------------------------------------------------------------------------------------------------------------------------
# Merge sortedName.bam file from each lane of the sample in one single file
#----------------------------------------------------------------------------------------------------------------------------
NUM_BAMS=$(ls -f $WORKDIR/BamMerged/ | wc -l);

if [ "$NUM_BAMS" -eq 1 ]; then
	cp $WORKDIR/BamMerged/${ID_SAMPLE}_${LANE}_to_${virus_name}_sortedName.bam $WORKDIR/${SAMPLE}_to_${virus_name}_allLanes.bam
else
	cd $WORKDIR/BamMerged
	samtools merge -@ $TH $WORKDIR/${SAMPLE}_to_${virus_name}_allLanes.bam *.bam
fi

#----------------------------------------------------------------------------------------------------------------------------
# Sorting resulting bam file and indexing for IGV visualization
#----------------------------------------------------------------------------------------------------------------------------
samtools sort -@ $TH -o $WORKDIR/${SAMPLE}_to_${virus_name}_sortedIGV.bam $WORKDIR/${SAMPLE}_to_${virus_name}_allLanes.bam
samtools index $WORKDIR/${SAMPLE}_to_${virus_name}_sortedIGV.bam

#----------------------------------------------------------------------------------------------------------------------------
# Sorting resulting bam file for assembly and Extract reads 1 and reads 2 
#----------------------------------------------------------------------------------------------------------------------------
samtools sort -@ $TH -n -o $WORKDIR/${SAMPLE}_to_${virus_name}_sortedName.bam $WORKDIR/${SAMPLE}_to_${virus_name}_allLanes.bam
samtools view -@ $TH -h -o $WORKDIR/${SAMPLE}_to_${virus_name}_sortedName.sam $WORKDIR/${SAMPLE}_to_${virus_name}_sortedName.bam

bedtools bamtofastq -i $WORKDIR/${SAMPLE}_to_${virus_name}_sortedName.bam \
-fq $WORKDIR/Fastq/${SAMPLE}_to_${virus_name}_Reads1.fq \
-fq2 $WORKDIR/Fastq/${SAMPLE}_to_${virus_name}_Reads2.fq

#----------------------------------------------------------------------------------------------------------------------------
# Assembly
#----------------------------------------------------------------------------------------------------------------------------
$TRINITY_EXE \
--CPU $TH --max_memory $MEM \
--seqType fq \
--left $WORKDIR/Fastq/${SAMPLE}_to_${virus_name}_Reads1.fq \
--right $WORKDIR/Fastq/${SAMPLE}_to_${virus_name}_Reads2.fq \
--output $WORKDIR/trinity \
--full_cleanup

mv $WORKDIR/trinity.Trinity.fasta $WORKDIR/${SAMPLE}_to_${virus_name}_assembly.fa
mv $WORKDIR/trinity.Trinity.fasta.gene_trans_map $WORKDIR/${SAMPLE}_to_${virus_name}_assembly_gene_trans.fa

#----------------------------------------------------------------------------------------------------------------------------
# Remove unused files
#----------------------------------------------------------------------------------------------------------------------------	
rm $WORKDIR/trinity

echo -e 'Time required'
END=$(date +%s);
echo $((END-START)) | awk '{print int($1/60)":"int($1%60)" min:sec"}'
echo -e 'End'
