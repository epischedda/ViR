import argparse
    
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
                                
def print_output_files(ReadsDict):
    '''Print output files'''
    print('Print intermediate files...')

    virus=[]
    virus_list=''
    reads=''
    readsID_searchInSam=''
    
    out_selected_reads_file_name=str(opts.outputPath)+'/Selected_Reads/Selected_Reads_Info.txt'
    out_selected_reads=open(out_selected_reads_file_name,'w')
    out_selected_reads.write('# SAMPLE_ID\tREAD_ID\tHR_CHR\tHR_START\tHR_MQ\tVR_SEQ\tVIRUS_ID\tVIRUS_START\tVIRUS_END\tVIRUS_SEQ\n')
                   
    for read_id, pair in ReadsDict.iteritems():
        SAMPLE_ID=str(pair.sample_name)
        READ_ID=str(pair.read_id)
        HR_CHR=str(pair.host_mate.chromosome)
        HR_START=str(pair.host_mate.start)
        HR_MQ=str(pair.host_mate.MQ)
        VR_SEQ=str(pair.viral_mate.viral_read_seq)
        VIRUS_ID=str(pair.viral_mate.virus_id)
        VIRUS_START=str(pair.viral_mate.virus_start)
        VIRUS_END=str(pair.viral_mate.virus_end)
        VIRAL_SEQ=str(pair.viral_mate.virus_seq)
  
        txt_line=SAMPLE_ID+'\t'+READ_ID+'\t'+HR_CHR+'\t'+HR_START+'\t'+HR_MQ+'\t'+VR_SEQ+'\t'+VIRUS_ID+'\t'+VIRUS_START+'\t'+VIRUS_END+'\t'+VIRAL_SEQ+'\n'        
        out_selected_reads.write(txt_line)

        readsID_searchInSam=readsID_searchInSam+str(READ_ID)+'\|'
        
    out_selected_reads.close()
    
    output_IDreads_listSam_file_name=str(opts.outputPath)+'/Selected_Reads/IDreads_listSam.txt'
    out_file_IDreads_listSam=open(output_IDreads_listSam_file_name,'w')
    
    if readsID_searchInSam is not '':
        readsID_searchInSam_list=readsID_searchInSam.split('\|')
        while len(readsID_searchInSam_list) > 500:
            readsIDtoprint='\|'.join(readsID_searchInSam_list[0:500])+'\n'
            readsID_searchInSam_list=readsID_searchInSam_list[500:len(readsID_searchInSam_list)]
            out_file_IDreads_listSam.write(readsIDtoprint)

        readsIDtoprint='\|'.join(readsID_searchInSam_list[0:-1])+'\n'
        out_file_IDreads_listSam.write(readsIDtoprint)
        
    out_file_IDreads_listSam.close()


def read_file_reader(txt):
    '''Read chimeric reads file'''
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
            
            VIRAL_SEQ=str(viral_read.virus_seq)
            lenViralSeq=len(VIRAL_SEQ)
            numA=VIRAL_SEQ.count('A')
            numT=VIRAL_SEQ.count('T')
            numC=VIRAL_SEQ.count('C')
            numG=VIRAL_SEQ.count('G')
            percAT=float(int(numA)+int(numT))/float(lenViralSeq)
            percAC=float(int(numA)+int(numC))/float(lenViralSeq)
            percAG=float(int(numA)+int(numG))/float(lenViralSeq)
            percCG=float(int(numG)+int(numC))/float(lenViralSeq)
            percCT=float(int(numT)+int(numC))/float(lenViralSeq)
            percGT=float(int(numT)+int(numG))/float(lenViralSeq)
            
            if (percAT) >= float(opts.percentage) or (percAC) >= float(opts.percentage) or (percAG) >= float(opts.percentage) or (percCG) >= float(opts.percentage) or (percCT) >= float(opts.percentage) or (percGT) >= float(opts.percentage) or int(lenViralSeq) <= int(opts.min_viral_length):
                continue
            else:
                ReadsDict[pair.read_id]=pair

    return ReadsDict


def main():
    parser = argparse.ArgumentParser('Filter couples of reads based on the viral sequence length and percentage of dinucleotides')
    parser.add_argument('-i_reads', '--file_reads', help="chimeric read file")
    parser.add_argument('-o', '--outputPath', help="output directory")
    parser.add_argument('-p', '--percentage', help="maximum percentage of each couple of dinucleotides accepted in a viral sequence", default=0.8)
    parser.add_argument('-min_vl', '--min_viral_length', help="minimum viral sequence length required to maintain a reads pair", default=30)
    
    global opts
    opts = parser.parse_args()
    
    read_file = open(opts.file_reads)
    ReadsDict=read_file_reader(read_file)
    
    print_output_files(ReadsDict)
    
main()