import argparse

class Pair_Info():
    sample_name=''
    read_id=''
    host_mate=[]
    viral_mate=[]
    
class HostRead_Info():
    chromosome=''
    start=0
    end=0
    host_read_seq=''
    flag=''
    ntBQ20=0
    MQ=0
    
class ViralRead_Info():
    viral_read_seq=''
    flag=''
    ntBQ20=0
    virus_id=''
    virus_start=0
    virus_end=0
    virus_seq=''
    virus_seq_len=0
    alignToRef='No'

def print_information(read_dict):
    '''print the output file'''  
    print('Print output files...')
    
    output_file_ChimericPairs_reads=str(opts.outputPath)+'/Final_ChimericPairs_Info.txt'
    out_file_ChimericPairs_reads=open(output_file_ChimericPairs_reads,'w')
    out_file_ChimericPairs_reads.write('# SAMPLE_ID\tREAD_ID\tHR_CHR\tHR_START\tHR_END\tHR_SEQ\tHR_FLAG\tHR_NT_BQ20\tHR_MQ\tVR_SEQ\tVR_FLAG\tVR_NT_BQ20\tVIRUS_ID\tVIRUS_START\tVIRUS_END\tVIRUS_SEQ\tVIRUS_SEQ_LEN\tVR_AlignToRef?\n')
    
    for read_id, pair in read_dict.iteritems():

        SAMPLE_ID=str(pair.sample_name)
        READ_ID=str(pair.read_id)

        HR_CHR=str(pair.host_mate.chromosome)
        HR_START=str(pair.host_mate.start)
        HR_END=str(pair.host_mate.end)
        HR_SEQ=str(pair.host_mate.host_read_seq)
        HR_FLAG=str(pair.host_mate.flag)
        HR_NT_BQ20=str(pair.host_mate.ntBQ20)
        HR_MQ=str(pair.host_mate.MQ)

        VR_SEQ=str(pair.viral_mate.viral_read_seq)
        VR_FLAG=str(pair.viral_mate.flag)
        VR_NT_BQ20=str(pair.viral_mate.ntBQ20)
        VIRUS_ID=str(pair.viral_mate.virus_id)
        VIRUS_START=str(pair.viral_mate.virus_start)
        VIRUS_END=str(pair.viral_mate.virus_end)
        VIRAL_SEQ=str(pair.viral_mate.virus_seq)
        VIRUS_SEQ_LEN=str(pair.viral_mate.virus_seq_len)
        VR_AlignToRef=str(pair.viral_mate.alignToRef)

        txt_line=SAMPLE_ID+'\t'+READ_ID+'\t'+HR_CHR+'\t'+HR_START+'\t'+HR_END+'\t'+HR_SEQ+'\t'+HR_FLAG+'\t'+HR_NT_BQ20+'\t'+HR_MQ+'\t'+VR_SEQ+'\t'+VR_FLAG+'\t'+VR_NT_BQ20+'\t'+VIRUS_ID+'\t'+VIRUS_START+'\t'+VIRUS_END+'\t'+VIRAL_SEQ+'\t'+VIRUS_SEQ_LEN+'\t'+VR_AlignToRef+'\n'
        out_file_ChimericPairs_reads.write(txt_line)

    out_file_ChimericPairs_reads.close()

        
def read_blastVR(ReadsDict, vr2ref_file):
    '''Read the blastn output file of the viral mate alignment to reference genome'''
    
    for line in vr2ref_file:

        line=line.rstrip()
        row=line.split('\t')
        
        ReadID=row[0]
        if ReadID in ReadsDict.keys():
            ReadsDict[ReadID].viral_mate.alignToRef='Yes'
        
    return ReadsDict

        
def read_uncomplete_result(txt):
    '''Read file from RC_MergeInfo.py'''
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
            host_read.end=row[4]
            host_read.host_read_seq=row[5]
            host_read.flag=row[6]
            host_read.ntBQ20=row[7]
            host_read.MQ=row[8]

            viral_read.viral_read_seq=row[9]
            viral_read.flag=row[10]
            viral_read.ntBQ20=row[11]
            viral_read.virus_id=row[12]
            viral_read.virus_start=row[13]
            viral_read.virus_end=row[14]
            viral_read.virus_seq=row[15]
            viral_read.virus_seq_len=row[16]

            pair.host_mate=host_read
            pair.viral_mate=viral_read
            
            ReadsDict[pair.read_id]=pair
    
    return ReadsDict
    

def main():
    parser = argparse.ArgumentParser('Generate the output files ViR_RefineCandidates.py')
    parser.add_argument('-i_vr2ref', '--file_vr2ref', help="blastn outfmt 6 output file")
    parser.add_argument('-i_chimeric_pair_info', '--file_chimeric_pair_info', help="file of the chimeric reads information")
    parser.add_argument('-o', '--outputPath', help="output directory")
    
    global opts
    opts = parser.parse_args()
    
    uncomplete_result_file = open(opts.file_chimeric_pair_info)
    ReadsDict=read_uncomplete_result(uncomplete_result_file)
    
    vr2ref_file = open(opts.file_vr2ref)
    ReadsDict=read_blastVR(ReadsDict, vr2ref_file)

    print_information(ReadsDict)
    
main()