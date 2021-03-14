#!/bin/bash -x

#----------------------------------------------------------------------------------------------------------------------------
# Input settings
#----------------------------------------------------------------------------------------------------------------------------
while [[ $# > 0 ]]
do
key="$1"

case $key in
	-work_files_dir)
	WORKDIR="$2"
	shift # past argument
	;;
	-sample_name)
	SAMPLE_NAME="$2"
	shift # past argument
	;;
	-sam_file)
	SAM_FILE="$2"
	shift # past argument
	;;
	-chimeric_reads_file)
	CHIMERIC_READS_FILE="$2"
	shift # past argument
	;;
	-out)
	OUTPUTDIR="$2"
	shift # past argument
	;;
	-max_percentage_dinucleotide_in_ViralSeq)
	MAX_PERC="$2"
	shift # past argument
	;;
	-minimum_virus_len)
	MIN_VIRUS_LEN="$2"
	shift # past argument
	;;
	-blastn_evalue)
	EVALUE="$2"
	shift # past argument
	;;
	-min_mate_distance)
	MATE_DISTANCE="$2"
	shift # past argument
	;;
	-reference_fasta)
	REF_FASTA="$2"
	shift # past argument
	;;
	-path_to_blastn)
	PATH_TO_BLASTN="$2"
	shift # past argument
	;;
	-path_to_bedtools)
	PATH_TO_BEDTOOLS="$2"
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
if [ -z "$MAX_PERC" ];then
  echo "Maximum percentage of dinucleotide in the viral sequence accepted unset, default value: 0.8"
  MAX_PERC=0.8
fi
if [ -z "$MIN_VIRUS_LEN" ];then
  echo "Minimum viral sequence lenght accepted unset, default value: 30"
  MIN_VIRUS_LEN=30
fi
if [ -z "$EVALUE" ];then
  echo "Maximum evalue accepted in blastn alignment of host and viral reads in the reference genome unset, default value: 1e-15"
  EVALUE=1e-15
fi
if [ -z "$MATE_DISTANCE" ];then
  echo "Minimum host-viral reads coordinates distance to mantain the pair as candidate unset, default value: 10000"
  MATE_DISTANCE=10000
fi

#----------------------------------------------------------------------------------------------------------------------------
# Creation of the directories
#----------------------------------------------------------------------------------------------------------------------------
if [ ! -d $OUTPUTDIR/$SAMPLE_NAME ]; then
	mkdir $OUTPUTDIR/$SAMPLE_NAME
	mkdir $OUTPUTDIR/$SAMPLE_NAME/OutputFiles
	mkdir $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles
	mkdir $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads
	printf "\n"
	echo -e $OUTPUTDIR/$SAMPLE_NAME 'Created'
fi

#----------------------------------------------------------------------------------------------------------------------------
# Select the header from the alignment
#----------------------------------------------------------------------------------------------------------------------------
python $WORKDIR/RC_HeaderSam.py -sam $SAM_FILE > $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/Selected_Reads.sam

#----------------------------------------------------------------------------------------------------------------------------
# Reads are filtered based on the viral sequence detected in the viral reads. In particular, viral sequences with more then 
# MAX_PERC of each dinucleotide combination and length less then MIN_VIRUS_LEN are filtered.
#----------------------------------------------------------------------------------------------------------------------------
echo -e 'Reads filtering...'
python $WORKDIR/RC_ReadsFiltering.py \
	-i_reads $CHIMERIC_READS_FILE \
	-o $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles \
	-p $MAX_PERC \
	-min_vl $MIN_VIRUS_LEN

num_viral_reads=$(wc -l < $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/IDreads_listSam.txt)

if [ "$num_viral_reads" = "0" ]; then
	echo -e 'No viral reads pass filters. Pipeline STOP!'
	
else
	#----------------------------------------------------------------------------------------------------------------------------
	# SAM file rows related to the READS ID of the remaining reads are selected.
	#----------------------------------------------------------------------------------------------------------------------------
	echo -e 'Look for read alignments in the sam file...'
	mkdir $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen
	mkdir $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen
	mkdir $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest
	mkdir $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Final_Reads

	while IFS='' read -r line || [[ -n "$line" ]]; do

		grep -w $line $SAM_FILE >> $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/Selected_Reads.sam
			
	done < $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/IDreads_listSam.txt

	#----------------------------------------------------------------------------------------------------------------------------
	# Merge information of host and viral reads; creates their fasta files
	#----------------------------------------------------------------------------------------------------------------------------
	python $WORKDIR/RC_MergeInfo.py \
	-i_selreads_sam $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/Selected_Reads.sam \
	-i_selreads_info $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/Selected_Reads_Info.txt \
	-o $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles

	#----------------------------------------------------------------------------------------------------------------------------
	# blastn viral reads against the reference genome
	#----------------------------------------------------------------------------------------------------------------------------
	echo -e 'Blastn Viral Reads to the reference genome...'

	$PATH_TO_BLASTN \
	-db $REF_FASTA \
	-query $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/All_ViralReads.fasta \
	-outfmt 6 \
	-out $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef.txt \
	-evalue $EVALUE
	
	python $WORKDIR/RC_Blast2Bed.py \
	-i $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef.txt \
	-o $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef.bed

	$PATH_TO_BEDTOOLS sort -i $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef.bed > $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef_sorted.bed

	rm $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef.bed
	
	#----------------------------------------------------------------------------------------------------------------------------
	# blastn host reads against the reference genome
	#----------------------------------------------------------------------------------------------------------------------------
	echo -e 'Blastn Host Reads to the reference genome...'

	$PATH_TO_BLASTN \
	-db $REF_FASTA \
	-query $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads.fasta \
	-outfmt 6 \
	-out $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef.txt \
	-evalue $EVALUE

	python $WORKDIR/RC_Blast2Bed.py \
	-i $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef.txt \
	-o $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef.bed

	$PATH_TO_BEDTOOLS sort -i $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef.bed > $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef_sorted.bed

	rm $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef.bed

	#----------------------------------------------------------------------------------------------------------------------------
	# Exclude mates that can map in the reference genome inside the window decided by the user
	#----------------------------------------------------------------------------------------------------------------------------
	grep '>' $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads.fasta > $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/List_ReadsFromSam.txt
	sed -i -e 's/>//g' $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/List_ReadsFromSam.txt

	touch $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest/Mates_inWindow.txt

	while IFS='' read -r line || [[ -n "$line" ]]; do
		
		grep $line $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef_sorted.bed > $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef_sorted_specificReads.bed
		grep $line $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef_sorted.bed > $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef_sorted_specificReads.bed

		if [[ -s $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef_sorted_specificReads.bed ]];
			then
				$PATH_TO_BEDTOOLS closest -d \
				-a $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef_sorted_specificReads.bed \
				-b $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef_sorted_specificReads.bed > $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest/ClosestMates_specificReads.txt

				awk -v d="$MATE_DISTANCE" '$13 >= 0 && $13 < d' $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest/ClosestMates_specificReads.txt >> $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest/Mates_inWindow.txt
		fi

	done < $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/List_ReadsFromSam.txt

	if [[ -s $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest/ClosestMates_specificReads.txt ]];
		then
		rm $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest/ClosestMates_specificReads.txt
	fi

	if [[ -s $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef_sorted_specificReads.bed ]];
		then
		rm $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef_sorted_specificReads.bed
	fi

	if [[ -s $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef_sorted_specificReads.bed ]];
		then
		rm $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef_sorted_specificReads.bed
	fi
	
	echo -e 'Exclude mates that map inside the window decided by the user...'
	python $WORKDIR/RC_FilterMatesInWindow.py \
	-i_mates_in_window $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest/Mates_inWindow.txt \
	-i_readsFromSam $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/List_ReadsFromSam.txt \
	-o $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest

	while IFS='' read -r line || [[ -n "$line" ]]; do

		grep -w $line $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Selected_Reads/ChimericPairs_Info.txt >> $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Final_Reads/ChimericPairs_Info.txt
		grep -w $line $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef.txt >> $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Final_Reads/blast_ViralReads_ToRef.txt
		grep -w $line $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef.txt >> $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Final_Reads/blast_HostReads_ToRef.txt
		grep $line $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/blast_ViralReads_ToRef_sorted.bed >> $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Final_Reads/blast_ViralReads_ToRef_sorted.bed
		grep $line $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/blast_HostReads_ToRef_sorted.bed >> $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Final_Reads/blast_HostReads_ToRef_sorted.bed
		grep -w $line -A 1 $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads.fasta >> $OUTPUTDIR/$SAMPLE_NAME/OutputFiles/Final_HostReads.fasta
		grep -w $line -A 1 $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Blastn_ViralReads_refgen/All_ViralReads.fasta >> $OUTPUTDIR/$SAMPLE_NAME/OutputFiles/Final_ViralReads.fasta
		
	done < $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Bedtools_closest/List_FinalReads.txt

	sed -i '/--/d' $OUTPUTDIR/$SAMPLE_NAME/OutputFiles/Final_HostReads.fasta
	sed -i '/--/d' $OUTPUTDIR/$SAMPLE_NAME/OutputFiles/Final_ViralReads.fasta

	#----------------------------------------------------------------------------------------------------------------------------
	# Generate output file
	#----------------------------------------------------------------------------------------------------------------------------
	python $WORKDIR/RC_GenerateOutput.py \
	-i_vr2ref $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Final_Reads/blast_ViralReads_ToRef.txt \
	-i_chimeric_pair_info $OUTPUTDIR/$SAMPLE_NAME/IntermediateFiles/Final_Reads/ChimericPairs_Info.txt \
	-o $OUTPUTDIR/$SAMPLE_NAME/OutputFiles
	
fi

printf "\n"     		
echo -e 'Time required'
END=$(date +%s);
echo $((END-START)) | awk '{print int($1/60)":"int($1%60)" min:sec"}'
echo -e 'End'
