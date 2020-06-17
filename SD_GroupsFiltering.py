import argparse
from collections import Counter
import os
from itertools import combinations
import math

class COMPLETE_INFO:
    sample_id=''
    read_id=''
    hr_chr=''
    hr_start=''
    hr_end=''
    hr_seq=''
    hr_flag=''
    hr_nt_bq20=''
    hr_mq=''
    vr_seq=''
    vr_flag=''
    vr_nt_bq20=''
    virus_id=''
    virus_start=''
    virus_end=''
    virus_seq=''
    virus_seq_len=0
    vr_alignToRef=''
    reads_group=[]
    second_alignment=''
    pass


class READS_GROUP:
    group_id=[]
    reads_list=[]
    second_align_chr=''
    second_align_start=''
    second_align_end=''
    ineve=''
    inPIRC=''
    nearTE=''
    virus_id=''
    virus_start=0
    virus_end=0
    max_virus_seq=0
    num_reads=0
    num_samples=0
    sec_alignment_subgroups={}
    num_equivalent_regions=0
    pass


def read_filebed(txt):
    '''Read file from filter bedtools closest'''
    print('Reads Region_Info')
    Coordinates_Dict={}

    for line in txt:
        if line.startswith('scaffold'):
            continue
        else:
            line=line.rstrip()
            parts=line.split('\t')
            
            scaffold=parts[0]
            start=parts[1]
            end=parts[2]
            id=parts[3]
            r_list=','.join(sorted(parts[5].split(',')))
            eve_name=parts[6]
            eve_start=parts[7]
            eve_end=parts[8]
            eve_distance=parts[9]
            piwic_name=parts[10]
            piwic_start=parts[11]
            piwic_end=parts[12]
            piwic_distance=parts[13]
            
            new_line=str(id)+'\t'+str(scaffold)+'\t'+str(start)+'\t'+str(end)+'\t'+str(eve_name)+'\t'+str(eve_start)+'\t'+str(eve_end)+'\t'+str(eve_distance)+'\t'+str(piwic_name)+'\t'+str(piwic_start)+'\t'+str(piwic_end)+'\t'+str(piwic_distance)+'\n'
            
            if r_list not in Coordinates_Dict.keys():
                hits_list=[]
                hits_list.append(new_line)
                Coordinates_Dict[r_list]=hits_list
            else:
                hits_list=Coordinates_Dict[r_list]
                hits_list.append(new_line)
                Coordinates_Dict[r_list]=hits_list

    print('Initial num read groups: '+str(len(Coordinates_Dict.keys())))

    return Coordinates_Dict


def read_filefasta(fasta):
    '''Read fasta file of host reads'''
    
    reads_list=[]

    for line in fasta:
        
        if line.startswith('>'):
            line=line.rstrip()
            parts=line.split('>')
            reads_list.append(str(parts[1]))
        
    return reads_list


def read_pair_info(in_all_pair_info):
    '''Read file of all chimeric reads'''
    
    Pair_Info_Dict={}

    for line in in_all_pair_info:
        
        if line.startswith('# SAMPLE_ID') or line.startswith('# STRAIN'):
		    continue
        else:
            line=line.rstrip()
            parts=line.split('\t')

            read_pair_info=COMPLETE_INFO()
            read_pair_info.sample_id=parts[0]
            read_pair_info.read_id=parts[1]
            read_pair_info.hr_chr=parts[2]
            read_pair_info.hr_start=parts[3]
            read_pair_info.hr_end=parts[4]
            read_pair_info.hr_seq=parts[5]
            read_pair_info.hr_flag=parts[6]
            read_pair_info.hr_nt_bq20=parts[7]
            read_pair_info.hr_mq=parts[8]
            read_pair_info.vr_seq=parts[9]
            read_pair_info.vr_flag=parts[10]
            read_pair_info.vr_nt_bq20=parts[11]
            read_pair_info.virus_id=parts[12]
            read_pair_info.virus_start=parts[13]
            read_pair_info.virus_end=parts[14]
            read_pair_info.virus_seq=parts[15]
            read_pair_info.virus_seq_len=parts[16]
            read_pair_info.vr_alignToRef=parts[17]
            Pair_Info_Dict[read_pair_info.read_id]=read_pair_info
           
    return Pair_Info_Dict


def read_group(Coordinates_Dict, reads_list, Pair_Info_Dict):
    '''1 - Groups in which reads have different virus id are splitted;
       2 - Control if some read groups is totally included in another, if yes consider the bigger read group; 
       3 - If two groups share at least a certain percentage (opts.perc_shared_reads) of the smaller of the 2 groups compared, they are merged and the positions of the bigger one are considered'''
    
    list_to_remove=[]
    db_to_insert={}

    # 1 - Groups in which reads have different virus_id are splitted
    print('\nGroups with different virus_id are splitted...')
    print('Num groups\tRemoved\tInserted')
    for k in Coordinates_Dict.iterkeys():
        id_reads=k.split(',')
        virus=[]
        
        for r in id_reads:
            virus.append(Pair_Info_Dict[r].virus_id)
         
        if int(len(set(virus))) != int(1):

            list_to_remove.append(k)
            v_dict=Counter(virus)

            for v_k, v_v in v_dict.iteritems():
                if int(v_v) > int(1):
                    
                    indexes=[i for i, e in enumerate(virus) if e == v_k]
                    sel_reads=','.join(sorted([id_reads[i] for i in indexes]))
                    
                    if sel_reads not in Coordinates_Dict.keys():
                        if sel_reads not in db_to_insert.keys():
                            db_to_insert[sel_reads]=Coordinates_Dict[k]
                        else:
                            temp_values=db_to_insert[sel_reads]+Coordinates_Dict[k]
                            db_to_insert[sel_reads]=temp_values

                    else:
                        temp_values=Coordinates_Dict[sel_reads]+Coordinates_Dict[k]
                        Coordinates_Dict[sel_reads]=temp_values


    print(str(len(Coordinates_Dict.keys()))+'\t'+str(len(list_to_remove))+'\t'+str(len(db_to_insert.keys())))                 
    for i in list_to_remove:
        del Coordinates_Dict[i]

    for k,v in db_to_insert.iteritems():
        Coordinates_Dict[k]=v
        
    
    list_to_remove=[]
    db_to_insert={}
    control=1
    counts=0
    num_removed_included=0
    num_removed_shared=0

    print('\nControl remaining groups...')
    print('Step\tNum groups\tRemoved_because_included\tCombinations\tRemoved_because_shared\tInserted')
    while control != 0:
        num_groups=len(Coordinates_Dict.keys())
        counts=counts+1
        control=0
        
        for k,j in combinations(Coordinates_Dict.keys(),2):
            k_reads=k.split(',')
            j_reads=j.split(',')
            
            # 2 - Control if some read groups is totally included in another, if yes consider the bigger read group;
            if len(list(set(k_reads) & set(j_reads)))== min(len(k_reads), len(j_reads)):
                control = 1
                if len(k_reads) <= len(j_reads):
                    if k not in list_to_remove:
                        list_to_remove.append(k)
                        num_removed_included=num_removed_included+1
                else:
                    if j not in list_to_remove:
                        list_to_remove.append(j)
                        num_removed_included=num_removed_included+1
        
        for i in list_to_remove:
            del Coordinates_Dict[i]

        list_to_remove=[]
        
        for k,j in combinations(Coordinates_Dict.keys(),2):
            k_reads=k.split(',')
            j_reads=j.split(',')

            share_num_reads=0
            if len(k_reads)>= len(j_reads):
                share_num_reads=math.floor(float(len(k_reads))*float(opts.perc_shared_reads))
            else:
                share_num_reads=math.floor(float(len(j_reads))*float(opts.perc_shared_reads))

            # 3 - If two groups share at least a certain percentage (opts.perc_shared_reads) of the smaller of the 2 groups compared, they are merged and the positions of the bigger one are considered;
            if int(len(list(set(k_reads) & set(j_reads)))) >= int(share_num_reads):
                
                control = 1
                new_group = ','.join(sorted(list(set(k_reads) | set(j_reads))))
                if len(k_reads) >= len(j_reads):
                    new_positions = Coordinates_Dict[k]
                else: 
                    new_positions = Coordinates_Dict[j]
                
                if new_group not in db_to_insert.keys():
                    db_to_insert[new_group]=new_positions
                        
                if k not in list_to_remove:
                    num_removed_shared=num_removed_shared+1
                    list_to_remove.append(k)
                if j not in list_to_remove:
                    num_removed_shared=num_removed_shared+1
                    list_to_remove.append(j)
        
        print(str(counts)+'\t'+str(num_groups)+'\t'+str(num_removed_included)+'\t'+str(len(list(combinations(Coordinates_Dict.keys(),2))))+'\t'+str(num_removed_shared)+'\t'+str(len(db_to_insert.keys())))
        for i in list_to_remove:
            del Coordinates_Dict[i]
        for k,v in db_to_insert.iteritems():
            Coordinates_Dict[k]=v
           
        list_to_remove=[]
        db_to_insert={}
        num_removed_included=0
        num_removed_shared=0

    RG_Dict={}
    i=1

    print('\nNumber of final groups: '+str(len(Coordinates_Dict.keys())))
    for k,v in Coordinates_Dict.iteritems():
        RG=READS_GROUP()
        RG.group_id='RG'+str(i)
        RG.reads_list=k.split(',')
        RG.num_reads=len(RG.reads_list)
        RG.second_align_chr=v[0].split('\t')[1]
        RG.second_align_start=v[0].split('\t')[2]
        RG.second_align_end=v[0].split('\t')[3]
        RG.num_equivalent_regions=len(v)
        
        temp_eve_list=[]
        temp_PIRC_list=[]
        for vi in v:
            temp_vi=vi.split('\t')
            if temp_vi[7] not in ['.','.\n']:
                if int(temp_vi[7]) >= int(0):
                    temp_eve_list.append(temp_vi[4])
            if temp_vi[11] not in ['.','.\n']:
                if int(temp_vi[11]) == int(0):
                    temp_PIRC_list.append(temp_vi[8])
        
        RG.ineve=','.join(list(set(temp_eve_list)))
        RG.inPIRC=','.join(list(set(temp_PIRC_list)))
        RG_Dict[RG.group_id]=RG
        i=i+1
    
    db_list_group={}
    temp_virus_id=[]
    temp_virus_start=[]
    temp_virus_end=[]
    temp_virus_seq_len=[]
    temp_num_samples=[]

    for k,v in RG_Dict.iteritems():
        for read in v.reads_list:
            if read not in db_list_group.keys():
                db_list_group[read]=[]
            db_list_group[read].append(k)

            temp_virus_id.append(Pair_Info_Dict[read].virus_id)
            temp_virus_start.append(int(Pair_Info_Dict[read].virus_start))
            temp_virus_end.append(int(Pair_Info_Dict[read].virus_end))
            temp_virus_seq_len.append(int(Pair_Info_Dict[read].virus_seq_len))
            temp_num_samples.append(Pair_Info_Dict[read].sample_id)

        v.virus_id=''.join(list(set(temp_virus_id)))
        v.virus_start=min(temp_virus_start)
        v.virus_end=max(temp_virus_end)
        v.max_virus_seq=max(temp_virus_seq_len)
        v.num_samples=len(list(set(temp_num_samples)))
        v.sec_alignment_subgroups=Counter(temp_num_samples)

        temp_virus_id=[]
        temp_virus_start=[]
        temp_virus_end=[]
        temp_virus_seq_len=[]
        temp_num_samples=[]
        
    return Coordinates_Dict , RG_Dict , db_list_group 


def write_output(Coordinates_Dict,Pair_Info_Dict, RG_Dict, db_list_group):
    '''Write output files:
    1 - Complete read pairs info with the read group selected
    2 - Bed file of 1 equivalent region per read group
    3 - Read groups table info
    4 - Complete information about equivalent regions in the genome for each group
    5 - Combination of sample id and read group to implement the second alignment'''

    out_complete=open(opts.outputPath+'/Complete_Dataset_Info.txt','w')
    out_groups_bed=open(opts.outputPath+'/Equivalent_region_per_Read_Group.bed','w')
    out_groups_info=open(opts.outputPath+'/Complete_Read_Groups_Info.txt','w')
    out_equivalent_regions=open(opts.outputPath+'/All_Equivalent_Regions_per_Read_Group.txt','w')

    # 1 - Complete read pairs info with read group selected
    out_complete.write('# SAMPLE_ID\tREAD_ID\tHR_CHR\tHR_START\tHR_END\tHR_SEQ\tHR_FLAG\tHR_NT_BQ20\tHR_MQ\tVR_SEQ\tVR_FLAG\tVR_NT_BQ20\tVIRUS_ID\tVIRUS_START\tVIRUS_END\tVIRUS_SEQ\tVIRUS_SEQ_LEN\tVR_AlignToRef?\tRead_Group\n')
    sample_list=[]
    list_read_ungrouped=[]
    num_shared_reads=0
    for readkey in sorted(Pair_Info_Dict.iterkeys()):
        readValue=Pair_Info_Dict[readkey]
        SAMPLE_ID=str(readValue.sample_id)
        READ_ID=str(readValue.read_id)

        HR_CHR=str(readValue.hr_chr)
        HR_START=str(readValue.hr_start)
        HR_END=str(readValue.hr_end)
        HR_SEQ=str(readValue.hr_seq)
        HR_FLAG=str(readValue.hr_flag)
        HR_NT_BQ20=str(readValue.hr_nt_bq20)
        HR_MQ=str(readValue.hr_mq)

        VR_SEQ=str(readValue.vr_seq)
        VR_FLAG=str(readValue.vr_flag)
        VR_NT_BQ20=str(readValue.vr_nt_bq20)
        VIRUS_ID=str(readValue.virus_id)
        VIRUS_START=str(readValue.virus_start)
        VIRUS_END=str(readValue.virus_end)
        VIRAL_SEQ=str(readValue.virus_seq)
        VIRUS_SEQ_LEN=str(readValue.virus_seq_len)
        VR_AlignToRef=str(readValue.vr_alignToRef)

        try:
            READ_GROUP=str(','.join(db_list_group[READ_ID]))
            readValue.reads_group=READ_GROUP
            if len(db_list_group[READ_ID]) > 1:
                num_shared_reads=num_shared_reads+1

        except:
            READ_GROUP='Ungrouped'
            readValue.reads_group=READ_GROUP
            list_read_ungrouped.append(READ_ID)
        
        txt_line=SAMPLE_ID+'\t'+READ_ID+'\t'+HR_CHR+'\t'+HR_START+'\t'+HR_END+'\t'+HR_SEQ+'\t'+HR_FLAG+'\t'+HR_NT_BQ20+'\t'+HR_MQ+'\t'+VR_SEQ+'\t'+VR_FLAG+'\t'+VR_NT_BQ20+'\t'+VIRUS_ID+'\t'+VIRUS_START+'\t'+VIRUS_END+'\t'+VIRAL_SEQ+'\t'+VIRUS_SEQ_LEN+'\t'+VR_AlignToRef+'\t'+READ_GROUP+'\n'
        
        if SAMPLE_ID not in sample_list:
            sample_list.append(SAMPLE_ID)
        out_complete.write(txt_line)

    out_complete.close()

    # 2 - Bed file of 1 equivalent region per read group
    for k in sorted(RG_Dict.iterkeys()):
        v=RG_Dict[k]
        out_groups_bed.write(str(v.second_align_chr)+'\t'+str(v.second_align_start)+'\t'+str(v.second_align_end)+'\t'+str(v.group_id)+'\n')
    
    out_groups_bed.close()

    # 3 - Read groups table info
    header_group_info='READ_GROUP\tSCAFFOLD\tSTART\tEND\tLENGHT\tNUM_EQ_RG\tnear/in_EVEs\tin_PIRCS\tnear_TE\tVIRUS_ID\tVIRUS_START\tVIRUS_END\tMAX_VIRAL_SEQ\tNUM_READS\tNUM_SAMPLES'
    for sample in sorted(sample_list):
        header_group_info=header_group_info+'\t'+str(sample)

    out_groups_info.write(header_group_info+'\n')

    for k in sorted(RG_Dict.iterkeys()):
        temp_value=RG_Dict[k]
        groups_info_line=''
        groups_info_line=str(temp_value.group_id)+'\t'+str(temp_value.second_align_chr)+'\t'+str(temp_value.second_align_start)+'\t'+str(temp_value.second_align_end)+'\t'+str(int(temp_value.second_align_end)-int(temp_value.second_align_start))+'\t'+str(int(temp_value.num_equivalent_regions)-1)+'\t'+str(temp_value.ineve)+'\t'+str(temp_value.inPIRC)+'\t'+str(temp_value.nearTE)+'\t'+str(temp_value.virus_id)+'\t'+str(temp_value.virus_start)+'\t'+str(temp_value.virus_end)+'\t'+str(temp_value.max_virus_seq)+'\t'+str(temp_value.num_reads)+'\t'+str(temp_value.num_samples)
        for sample in sorted(sample_list):
            try:
                groups_info_line=groups_info_line+'\t'+str(temp_value.sec_alignment_subgroups[sample])
            except:
                groups_info_line=groups_info_line+'\t'+str(0)

        out_groups_info.write(groups_info_line+'\n')
    
    virus_ungruped={}
    for r in sorted(list_read_ungrouped):
        read=Pair_Info_Dict[r]
        if read.virus_id not in virus_ungruped.keys():
            virus_ungruped[read.virus_id]=[read.sample_id]
        else:
            temp_sample=virus_ungruped[read.virus_id]
            temp_sample.append(read.sample_id)
            virus_ungruped[read.virus_id]=temp_sample

    for k in virus_ungruped.iterkeys():
        temp=Counter(virus_ungruped[k])
        groups_info_line=str('Ungrouped')+'\t\t\t\t\t\t\t\t\t'+str(k)+'\t\t\t\t'+str(len(virus_ungruped[k]))+'\t'+str(len(temp.keys()))
        
        for sample in sorted(sample_list):
            try:
                groups_info_line=groups_info_line+'\t'+str(temp[sample])
            except:
                groups_info_line=groups_info_line+'\t'+str(0)
        
        out_groups_info.write(groups_info_line+'\n')
    
    groups_info_line=str('Shared')+'\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(num_shared_reads)+'\t'
    out_groups_info.write(groups_info_line+'\n')    
    out_groups_info.close()

    # 4 - Complete information about equivalent regions in the genome for each group
    for k,v in Coordinates_Dict.iteritems():
        rgid=''
        for z in RG_Dict.itervalues():
            if ','.join(z.reads_list) == k:
                rgid=z.group_id
            
        main_line=str(rgid)+'\n'+str(len(k.split(',')))+' reads\n'+str(k.split(','))+'\n\n'
        def_line='id\tscaffold\tstart\tend\tEVE_name\tEVE_start\tEVE_end\tEVE_distance\tpiwi_cluster_name\tpiwi_cluster_start\tpiwi_cluster_end\tpiwi_cluster_distance\n'
        out_equivalent_regions.write(main_line+def_line)
        
        for hit in v:
            out_equivalent_regions.write(hit)

        out_equivalent_regions.write('\n\n')

    out_equivalent_regions.close()

    # 5 - Combination of sample id and read group to implement the second alignment
    os.mkdir(opts.outputPath+'/RG/')
    for rg,v in RG_Dict.iteritems():
        if v.ineve=='':
            for sample,count_num in v.sec_alignment_subgroups.items():
                if int(count_num) > int(1):
                    if not os.path.exists(opts.outputPath+'/RG/'+str(rg)):
                        os.mkdir(opts.outputPath+'/RG/'+str(rg))
                    if not os.path.exists(opts.outputPath+'/RG/'+str(rg)+'/'+str(sample)):
                        os.mkdir(opts.outputPath+'/RG/'+str(rg)+'/'+str(sample))
                    
                    out_rlist=open(opts.outputPath+'/RG/'+str(rg)+'/'+str(sample)+'/Reads.list','w')
                    for read in sorted(Pair_Info_Dict.iterkeys()):
                        readValue=Pair_Info_Dict[read]
                        if str(sample) == str(readValue.sample_id):
                            for r in readValue.reads_group.split(','):
                                if str(rg) == str(r):
                                    out_rlist.write(str(readValue.read_id)+'\n')
                    out_rlist.close()

                    print(rg+'\t'+sample+'\t'+str(count_num))

    
def main():
    parser = argparse.ArgumentParser('Select read groups for each read id; select for each group all the possible positions in the reference genome in which the group of reads can map (equivalent regions)')
    parser.add_argument('-i', '--filebed', help="file from bedtools closest")
    parser.add_argument('-i_reads', '--all_host_fasta', help="fasta file of the host reads to extract unique reads list")
    parser.add_argument('-i_pair_info', '--all_pair_info_txt', help="Chimeric Pair Information from ViR_RefineCandidates")
    parser.add_argument('-perc_shared_reads', '--perc_shared_reads', help="minimum percentage of reads that two groups have to share to be merged", default=0.8)
    parser.add_argument('-o', '--outputPath', help="absolute path for output files")

    global opts
    opts = parser.parse_args()
    
    in_all_pair_info = open(opts.all_pair_info_txt)
    Pair_Info_Dict = read_pair_info(in_all_pair_info)
    
    in_file = open(opts.filebed)
    Coordinates_Dict=read_filebed(in_file)
    
    in_fasta = open(opts.all_host_fasta)
    reads_list=read_filefasta(in_fasta)
    
    [Coordinates_Dict,RG_Dict,db_list_group]=read_group(Coordinates_Dict,reads_list,Pair_Info_Dict)
    write_output(Coordinates_Dict,Pair_Info_Dict,RG_Dict,db_list_group)
    
main()