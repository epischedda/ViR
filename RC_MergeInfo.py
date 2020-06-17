import argparse
    
class Mate():
    ReadID=''
    Flag=''
    Contig=''
    PosStart=0
    MQ=0
    CIGAR=''
    RefNAMEother_mate=''
    PosNAMEother_mate=0
    TLEN=0
    SEQ=''
    lenNtQUAL20=''

class Pair_Sam():
    Host_Read=Mate()
    Viral_Read=Mate()

class Pair_Info():
    sample_name=''
    read_id=''
    host_mate=[]
    viral_mate=[]
    
class HostRead_Info():
    chromosome=''
    start=0
    MQ=0
    
class ViralRead_Info():
    viral_read_seq=''
    virus_id=''
    virus_start=0
    virus_end=0
    virus_seq=''


def print_information(read_dict, align_dict):
    '''print output files'''  
    
    output_file_ChimericPairs_reads=str(opts.outputPath)+'/Selected_Reads/ChimericPairs_Info.txt'
    out_file_ChimericPairs_reads=open(output_file_ChimericPairs_reads,'w')
    out_file_ChimericPairs_reads.write('# SAMPLE_ID\tREAD_ID\tHR_CHR\tHR_START\tHR_END\tHR_SEQ\tHR_FLAG\tHR_NT_BQ20\tHR_MQ\tVR_SEQ\tVR_FLAG\tVR_NT_BQ20\tVIRUS_ID\tVIRUS_START\tVIRUS_END\tVIRUS_SEQ\tVIRUS_SEQ_LEN\n')
    
    output_file_viral_fasta=str(opts.outputPath)+'/Blastn_ViralReads_refgen/All_ViralReads.fasta'
    out_file_viral_fasta=open(output_file_viral_fasta,'w')
    
    output_file_host_fasta=str(opts.outputPath)+'/Blastn_HostReads_refgen/All_HostReads.fasta'
    out_file_host_fasta=open(output_file_host_fasta,'w')
    
    for read_id, pair in read_dict.iteritems():
        val=align_dict[read_id]

        SAMPLE_ID=str(pair.sample_name)
        READ_ID=str(pair.read_id)

        HR_CHR=str(pair.host_mate.chromosome)
        HR_START=str(pair.host_mate.start)
        HR_SEQ=str(val.Host_Read.SEQ)
        HR_END=str(int(HR_START)+int(len(HR_SEQ)))
        HR_FLAG=str(val.Host_Read.Flag)
        HR_NT_BQ20=str(val.Host_Read.lenNtQUAL20)
        HR_MQ=str(pair.host_mate.MQ)

        VR_SEQ=str(pair.viral_mate.viral_read_seq)
        VR_FLAG=str(val.Viral_Read.Flag)
        VR_NT_BQ20=str(val.Viral_Read.lenNtQUAL20)
        VIRUS_ID=str(pair.viral_mate.virus_id)
        VIRUS_START=str(pair.viral_mate.virus_start)
        VIRUS_END=str(pair.viral_mate.virus_end)
        VIRAL_SEQ=str(pair.viral_mate.virus_seq)
        VIRUS_SEQ_LEN=str(len(VIRAL_SEQ))

        txt_line=SAMPLE_ID+'\t'+READ_ID+'\t'+HR_CHR+'\t'+HR_START+'\t'+HR_END+'\t'+HR_SEQ+'\t'+HR_FLAG+'\t'+HR_NT_BQ20+'\t'+HR_MQ+'\t'+VR_SEQ+'\t'+VR_FLAG+'\t'+VR_NT_BQ20+'\t'+VIRUS_ID+'\t'+VIRUS_START+'\t'+VIRUS_END+'\t'+VIRAL_SEQ+'\t'+VIRUS_SEQ_LEN+'\n'
        out_file_ChimericPairs_reads.write(txt_line)            
        

        string_header_fasta='>'+READ_ID+'\n'
        out_file_viral_fasta.write(string_header_fasta+VR_SEQ+'\n')
        out_file_host_fasta.write(string_header_fasta+HR_SEQ+'\n')
        
    out_file_ChimericPairs_reads.close()
    out_file_viral_fasta.close()
    out_file_host_fasta.close()
    

def read_alignment_sam(txt, ReadsDict):
    '''Read the sam file of the selected chimeric reads'''    
    AlignmentDict={}
    
    for r in ReadsDict.keys():
        pair=Pair_Sam()
        AlignmentDict[r]=pair
        
    for line in txt:
        if line.startswith('@'):
            continue
        else:
            viral_mate_sam=Mate()
            host_mate_sam=Mate()
            
            line=line.rstrip()
            row=line.split('\t')
            
            READID=row[0]
            SEQ=row[9]

            if READID in ReadsDict.keys():
                SEQ_VIRAL_READ=ReadsDict[READID].viral_mate.viral_read_seq
                
                if str(SEQ)==str(SEQ_VIRAL_READ):
                    viral_mate_sam.ReadID=READID
                    viral_mate_sam.Flag=row[1]
                    viral_mate_sam.Contig=row[2]
                    viral_mate_sam.PosStart=row[3]
                    viral_mate_sam.MQ=row[4]
                    viral_mate_sam.CIGAR=row[5]
                    viral_mate_sam.RefNAMEother_mate=row[6]
                    viral_mate_sam.PosNAMEother_mate=row[7]
                    viral_mate_sam.TLEN=row[8]
                    viral_mate_sam.SEQ=SEQ
                    viral_mate_sam.lenNtQUAL20=len(row[10])-(row[10].count('!')+row[10].count('"')+row[10].count('#')+row[10].count('$')+row[10].count('%')+row[10].count('&')+row[10].count('\'')+row[10].count('(')+row[10].count(')')+row[10].count('*')+row[10].count('+')+row[10].count(',')+row[10].count('-')+row[10].count('.')+row[10].count('/')+row[10].count('0')+row[10].count('1')+row[10].count('2')+row[10].count('3')+row[10].count('4')+row[10].count('5'))               
                    
                    al_dict_val=AlignmentDict[READID]
                    al_dict_val.Viral_Read=viral_mate_sam            
                    AlignmentDict[READID]=al_dict_val
                    
                else:
                    host_mate_sam.ReadID=READID
                    host_mate_sam.Flag=row[1]
                    host_mate_sam.Contig=row[2]
                    host_mate_sam.PosStart=row[3]
                    host_mate_sam.MQ=row[4]
                    host_mate_sam.CIGAR=row[5]
                    host_mate_sam.RefNAMEother_mate=row[6]
                    host_mate_sam.PosNAMEother_mate=row[7]
                    host_mate_sam.TLEN=row[8]
                    host_mate_sam.SEQ=SEQ
                    host_mate_sam.lenNtQUAL20=len(row[10])-(row[10].count('!')+row[10].count('"')+row[10].count('#')+row[10].count('$')+row[10].count('%')+row[10].count('&')+row[10].count('\'')+row[10].count('(')+row[10].count(')')+row[10].count('*')+row[10].count('+')+row[10].count(',')+row[10].count('-')+row[10].count('.')+row[10].count('/')+row[10].count('0')+row[10].count('1')+row[10].count('2')+row[10].count('3')+row[10].count('4')+row[10].count('5'))
                    
                    al_dict_val=AlignmentDict[READID]
                    al_dict_val.Host_Read=host_mate_sam            
                    AlignmentDict[READID]=al_dict_val
        
    return AlignmentDict;        
     
            
def read_complete_result(txt):
    '''Read the Selected_Reads_Info.txt file from RC_ReadsFiltering.py'''
    ReadsDict={}
    
    for line in txt:
        pair=Pair_Info()
        host_read=HostRead_Info()
        viral_read=ViralRead_Info()
    
        if line.startswith('#'):
            continue   
        else:
            
            line=line.rstrip()
            row=line.split('\t')
            
            pair.sample_name=row[0]
            pair.read_id=row[1]
            
            host_read.chromosome=row[2]
            host_read.start=row[3]
            host_read.MQ=row[4]

            viral_read.viral_read_seq=row[5]
            viral_read.virus_id=row[6]
            viral_read.virus_start=row[7]
            viral_read.virus_end=row[8]
            viral_read.virus_seq=row[9]

            pair.host_mate=host_read
            pair.viral_mate=viral_read
            
            ReadsDict[pair.read_id]=pair
    
    return ReadsDict
    

def main():
    parser = argparse.ArgumentParser('Merges information of host and viral reads; creates their fasta files.')
    parser.add_argument('-i_selreads_sam', '--file_sam', help="sam file of the selected chimeric reads")
    parser.add_argument('-i_selreads_info', '--file_info', help="text file of the selected viral reads")
    parser.add_argument('-o', '--outputPath', help="output directory")
    
    global opts
    opts = parser.parse_args()
    
    selected_reads_info = open(opts.file_info)
    ReadsDict=read_complete_result(selected_reads_info)

    sam_file = open(opts.file_sam)
    AlignmentDict=read_alignment_sam(sam_file, ReadsDict)

    print_information(ReadsDict, AlignmentDict)
    
main()