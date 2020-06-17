#!/bin/bash -x

#----------------------------------------------------------------------------------------------------------------------------
# Input settings
#----------------------------------------------------------------------------------------------------------------------------

while [[ $# > 0 ]]
do
key="$1"

case $key in
  -analysis_name)
  ANALYSIS_NAME="$2"
  shift # past argument
  ;;
  -sample_list)
  SAMPLE_LIST="$2"
  shift # past argument
  ;;
  -work_files_dir)
  WORKDIR="$2"
  shift # past argument
  ;;
  -outdict_RefCand)
  REFCAND_DIR="$2"
  shift # past argument
  ;;
  -out)
  OUTPUTDIR="$2"
  shift # past argument
  ;;
  -reference_fasta)
  REF_FASTA="$2"
  shift # past argument
  ;;
  -repreg_fasta)
  REPREG_FASTA="$2"
  shift # past argument
  ;;
  -bed_EVE_annotated)
  BED_EVE="$2"
  shift # past argument
  ;;
  -bed_piwi_clusters)
  PIWI_CLUSTERS="$2"
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
  -trinity_exe)
  PATH_TO_TRINITY="$2"
  shift # past argument
  ;;
  -samtools_exe)
  PATH_TO_SAMTOOLS="$2"
  shift # past argument
  ;;
  -bwa_exe)
  PATH_TO_BWA="$2"
  shift # past argument
  ;;
  -merge_dist)
  MERGE_DIST="$2"
  shift # past argument
  ;;
  -eve_dist)
  EVE_DIST="$2"
  shift # past argument
  ;;
  -piwi_dist)
  PIWI_DIST="$2"
  shift # past argument
  ;;
  -minReads_inRegion)
  MIN_READS="$2"
  shift # past argument
  ;;
  -percReadsShared_inGroup_union)
  PERC_READS_SHARED="$2"
  shift # past argument
  ;;
  -min_TE_al_length)
  MIN_AL_LEN="$2"
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
if [ -z "$PERC_READS_SHARED" ];then
  echo "Minimum percentage of reads shared among group to determine their union unset, default value: 0.8"
  PERC_READS_SHARED=0.8
fi
if [ -z "$MIN_READS" ];then
  echo "Minimum number of reads within a group unset, default value: 2"
  MIN_READS=2
fi
if [[ ( ! -z "$PIWI_CLUSTERS" ) && ( -z "$PIWI_DIST" )]];then
  echo "Maximum distance between candidates coordinates and piwi RNA clusters unset, default value: 1"
  PIWI_DIST=1
fi
if [[ ( ! -z "$BED_EVE" ) && ( -z "$EVE_DIST" )]];then
  echo "Maximum distance between candidates coordinates and EVEs annotated unset, default value: 10000"
  EVE_DIST=10000
fi
if [ -z "$MERGE_DIST" ];then
  echo "Maximum distance between reads coordinates to be considered as group unset, default value: 1000"
  MERGE_DIST=1000
fi
if [[ ( ! -z "$REPREG_FASTA" ) && ( -z "$MIN_AL_LEN" )]];then
  echo "Minimum alignment lenght accepted for TE in reference genome unset, default value: 100"
  MIN_AL_LEN=100
fi

#----------------------------------------------------------------------------------------------------------------------------
# Creation of the directories
#----------------------------------------------------------------------------------------------------------------------------

if [ ! -d $OUTPUTDIR/$ANALYSIS_NAME ]; then
	mkdir $OUTPUTDIR/$ANALYSIS_NAME
	mkdir $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles
	mkdir $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles
  printf "\n"
	echo -e $OUTPUTDIR/$ANALYSIS_NAME 'Created'
fi

#----------------------------------------------------------------------------------------------------------------------------
# Collect information from samples in sample.list file
#----------------------------------------------------------------------------------------------------------------------------
echo -e 'Collect information from samples...'
mkdir $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen
mkdir $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Merge_ChimericReadsPairInfo

while IFS='' read -r line || [[ -n "$line" ]]; do

  cat $REFCAND_DIR/$line/IntermediateFiles/Final_Reads/blast_HostReads_ToRef_sorted.bed >>  $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef.bed
  cat $REFCAND_DIR/$line/OutputFiles/Final_HostReads.fasta >>  $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads.fasta
  cat $REFCAND_DIR/$line/OutputFiles/Final_ChimericPairs_Info.txt >> $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Merge_ChimericReadsPairInfo/All_ChimericPairs_Info.txt

done < $SAMPLE_LIST

if [[ -s $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads.fasta ]] ; then
  
  #----------------------------------------------------------------------------------------------------------------------------
  # Merge mapping coordinates of the host reads in the reference genome, to select regions overlapped
  #----------------------------------------------------------------------------------------------------------------------------
  echo -e 'Merge mapping coordinates of the host reads in the reference genome...'

  $PATH_TO_BEDTOOLS sort -i $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef.bed \
  > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_sorted.bed

  rm $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef.bed

  $PATH_TO_BEDTOOLS merge -d $MERGE_DIST \
  -i $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_sorted.bed \
  -c 1,4 -o count,collapse \
  > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_sorted_merged.bed

  #----------------------------------------------------------------------------------------------------------------------------
  # Filter out merged coordinates in which host reads are less than $MIN_READS
  #----------------------------------------------------------------------------------------------------------------------------
  echo -e 'Filter out merged coordinates including few reads...'

  python $WORKDIR/SD_GroupsDefinition.py \
  -i $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_sorted_merged.bed \
  -o_bed $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_selectedCoordinates.bed \
  -o_region_info $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/RegionsInfo.txt \
  -minReads $MIN_READS

  if [ -z "$BED_EVE" ];then
    echo "Bed file of EVE annotated not available!"
  else
    #----------------------------------------------------------------------------------------------------------------------------
    # Control if there is some EVE annotated in the boundary regions ($EVE_DIST) of each result
    #----------------------------------------------------------------------------------------------------------------------------
    echo -e 'Control if already annotated EVEs are in the boundary regions...'

    mkdir $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated

    EVE_SORTED=$(echo $BED_EVE | rev | cut -f 1 -d'.' | rev);

    $PATH_TO_BEDTOOLS sort -i $BED_EVE > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/EVE_sorted.bed
    $PATH_TO_BEDTOOLS sort -i $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_selectedCoordinates.bed \
    > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_selectedCoordinates_sorted.bed

    $PATH_TO_BEDTOOLS closest -d \
    -a $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_selectedCoordinates_sorted.bed \
    -b $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/EVE_sorted.bed \
    > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/Regions_closestEVE.txt

    rm $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/EVE_sorted.bed
    rm $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_selectedCoordinates_sorted.bed
    
    if [[ -s $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/Regions_closestEVE.txt ]] ; then
      tr -s '\t' '\t' < $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/Regions_closestEVE.txt \
      > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/Regions_closestEVE_corrected.txt

      python $WORKDIR/SD_FilterBedtoolsClosest_byDistance.py \
      -i $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/Regions_closestEVE_corrected.txt \
      -o $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/Regions_closestEVE_corrected_filteredByDist.txt \
      -d $EVE_DIST \
      -seq_type EVE
    else
      BED_EVE=
    fi
  fi

  if [ -z "$PIWI_CLUSTERS" ];then
    echo "Bed file of piwi clusters not available!"
  else
    #----------------------------------------------------------------------------------------------------------------------------
    # Control if there is some PIWI_clusters annotated in the boundary regions ($MAX_DIST_CLOSEST_EVE) of each result
    #----------------------------------------------------------------------------------------------------------------------------
    echo -e 'Control if regions are in annotated piwi clusters...'

    mkdir $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters

    PIWICLUSTER_SORTED=$(echo $PIWI_CLUSTERS | rev | cut -f 1 -d'.' | rev);

    $PATH_TO_BEDTOOLS sort -i $PIWI_CLUSTERS > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/piwicluster_sorted.bed
    $PATH_TO_BEDTOOLS sort -i $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_selectedCoordinates.bed \
    > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_selectedCoordinates_sorted.bed

    $PATH_TO_BEDTOOLS closest -d \
    -a $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_selectedCoordinates_sorted.bed \
    -b $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/piwicluster_sorted.bed \
    > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/Regions_closestPiwiClusters.txt

    rm $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/piwicluster_sorted.bed
    rm $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads_ToRef_selectedCoordinates_sorted.bed

    if [[ -s $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/Regions_closestPiwiClusters.txt ]] ; then
      tr -s '\t' '\t' < $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/Regions_closestPiwiClusters.txt \
      > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/Regions_closestPiwiClusters_corrected.txt

      python $WORKDIR/SD_FilterBedtoolsClosest_byDistance.py \
      -i $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/Regions_closestPiwiClusters_corrected.txt \
      -o $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/Regions_closestPiwiClusters_corrected_filteredByDist.txt \
      -d $PIWI_DIST \
      -seq_type piwi

    else
      PIWI_CLUSTERS=
    fi
  fi

  #----------------------------------------------------------------------------------------------------------------------------
  # Join information from the boundary regions about EVE annatated and piwi clusters, if available
  #----------------------------------------------------------------------------------------------------------------------------

  if [[ ( ! -z "$PIWI_CLUSTERS" ) && ( ! -z "$BED_EVE" )]];then
      echo -e 'Add information of EVEs and piwi clusters...'

    join -1 4 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.7,2.8,2.9,2.10 \
    <(sort -k4 $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/Regions_closestEVE_corrected_filteredByDist.txt) \
    <(sort -k4 $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/Regions_closestPiwiClusters_corrected_filteredByDist.txt) \
    > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/RegionsInfo_final.txt
  fi
  if [[ -z "$PIWI_CLUSTERS" && ( ! -z "$BED_EVE" )]];then
    echo -e 'Add information of EVEs...'

    join -1 4 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.7,2.8,2.9,2.10 \
    <(sort -k4 $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_EVEannotated/Regions_closestEVE_corrected_filteredByDist.txt) \
    <(sort -k4 $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/RegionsInfo.txt) \
    > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/RegionsInfo_final.txt
  fi
  if [[ ( ! -z "$PIWI_CLUSTERS" ) && ( -z "$BED_EVE" )]];then
    echo -e 'Add information of piwi clusters...'

    join -1 4 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.7,2.8,2.9,2.10 \
    <(sort -k4 $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/RegionsInfo.txt) \
    <(sort -k4 $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Bedtools_closest_PiwiClusters/Regions_closestPiwiClusters_corrected_filteredByDist.txt) \
    > $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/RegionsInfo_final.txt

  fi
  if [[ ( -z "$PIWI_CLUSTERS" ) && ( -z "$BED_EVE" )]];then
    cp $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/RegionsInfo.txt $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/RegionsInfo_final.txt

  fi
  #----------------------------------------------------------------------------------------------------------------------------
  # Print final table
  #----------------------------------------------------------------------------------------------------------------------------
  echo -e 'Print final table...'

  python $WORKDIR/SD_GroupsFiltering.py \
  -i $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/RegionsInfo_final.txt \
  -i_reads $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_HostReads_refgen/All_HostReads.fasta \
  -i_pair_info $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Merge_ChimericReadsPairInfo/All_ChimericPairs_Info.txt \
  -perc_shared_reads $PERC_READS_SHARED \
  -o $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles

  #----------------------------------------------------------------------------------------------------------------------------
  # Second Alignment
  #----------------------------------------------------------------------------------------------------------------------------

  printf "\n"
  echo -e 'Second alignment...'
  $PATH_TO_BEDTOOLS sort -i $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Equivalent_region_per_Read_Group.bed \
  > $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Equivalent_region_per_Read_Group_sorted.bed

  rm $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Equivalent_region_per_Read_Group.bed

  $PATH_TO_BEDTOOLS getfasta -name \
  -fi $REF_FASTA \
  -bed $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Equivalent_region_per_Read_Group_sorted.bed \
  -fo $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Equivalent_region_per_Read_Group.fasta

  for gr in $( ls $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/RG ); do
    for smp in $( ls $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/RG/$gr ); do
      grep -A1 -w $gr $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Equivalent_region_per_Read_Group.fasta > $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/RG/$gr/$smp/$gr.fasta
      echo $gr $smp
    
      bash $WORKDIR/ViR_AlignToGroup.sh \
      -work_files_dir $WORKDIR \
      -sample_sam $REFCAND_DIR/$smp/IntermediateFiles/Selected_Reads/Selected_Reads.sam \
      -sample_name $smp \
      -group_name $gr \
      -groups_fasta $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/RG/$gr/$smp/$gr.fasta \
      -reads_list $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/RG/$gr/$smp/Reads.list \
      -outdir $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/RG/$gr/$smp \
      -bedtools_exe $PATH_TO_BEDTOOLS \
      -trinity_exe $PATH_TO_TRINITY \
      -samtools_exe $PATH_TO_SAMTOOLS \
      -bwa_exe $PATH_TO_BWA

      printf "\n" 
    done
  done

  #----------------------------------------------------------------------------------------------------------------------------
  # Blastn alignment of selected regions on the repeated region fasta file
  #----------------------------------------------------------------------------------------------------------------------------

  if [ -z "$REPREG_FASTA" ];then
    echo "Fasta file of repeated elements not available!"
  else
    
    mkdir $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_regions_repreg

    if [[ -s "$OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Equivalent_region_per_Read_Group.fasta" ]]; then
      $PATH_TO_BLASTN -db $REPREG_FASTA \
      -query $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Equivalent_region_per_Read_Group.fasta \
      -outfmt 6 \
      -out $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_regions_repreg/blastn_EqReg_onRepReg.txt \
      -task blastn

      python $WORKDIR/SD_AddTE_Info.py \
      --blastfile $OUTPUTDIR/$ANALYSIS_NAME/IntermediateFiles/Blastn_regions_repreg/blastn_EqReg_onRepReg.txt \
      --groupInfo $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Complete_Read_Groups_Info.txt \
      --min_len $MIN_AL_LEN \
      --outputPath $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles
      
      rm $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Complete_Read_Groups_Info.txt
      mv $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Complete_Read_Groups_Info_final.txt $OUTPUTDIR/$ANALYSIS_NAME/OutputFiles/Complete_Read_Groups_Info.txt
    else
      echo "No Read Groups detected!"
    fi  
  fi

else
  echo -e 'No reads to analyze'
fi

echo -e 'Time required'
END=$(date +%s);
echo $((END-START)) | awk '{print int($1/60)":"int($1%60)" min:sec"}'
echo -e 'End'