import argparse

def read_file_blast(txt):
    '''Read output table outfmt 6 of blastn and create a non sorted bed file: make unique sequence name in case of multiple alignment and put in order start<stop inserting strand column in bed file'''
    
    list_reads=[]
    blast_db={}
    line_name=''
    out_file_tot=open(opts.outputFile,'w')
    
    for line in txt:
        if line.startswith('#'):
            continue
        else:
            start=0
            end=0
            
            line=line.rstrip()
            parts=line.split('\t')
            query_acc_ver=parts[0]
            subject_acc_ver=parts[1]
            identity=parts[2]
            alignment_length=parts[3]
            s_start=parts[8]
            s_end=parts[9]
            
            if query_acc_ver not in list_reads:    
                n = 1
                list_reads.append(query_acc_ver)
            
            if str(query_acc_ver) != str(line_name):
                new_query_acc_ver = query_acc_ver+'_'+str(n)
                n = n+1
                line_name = new_query_acc_ver
                
            else:
                line_name=query_acc_ver
            
            if int(s_start) > int(s_end):
                start = str(s_end)
                end = str(s_start)
                strand = '-'
            else:
                start = str(s_start)
                end = str(s_end)
                strand = '+'

            if int(alignment_length) >= int(opts.minAlignLen) and float(identity) > float(opts.min_id) and float(identity) <= float(opts.max_id):
                out_line=str(subject_acc_ver)+'\t'+str(start)+'\t'+str(end)+'\t'+str(line_name)+'\t0\t'+strand+'\n'    
                out_file_tot.write(out_line)
    
    out_file_tot.close()
    return blast_db


def main():
    parser = argparse.ArgumentParser('Convert blast result table in bed format. Outputh bed file have to be sorted.')
    parser.add_argument('-i', '--file', help="blast result table as input file")
    parser.add_argument('-o', '--outputFile', help="output path file in bed format")
    parser.add_argument('-len', '--minAlignLen', help="minimum hit alignment length", default=0)
    parser.add_argument('-min_id', '--min_id', help="minimum hit identity accepted", default=0)
    parser.add_argument('-max_id', '--max_id', help="maximum hit identity accepted", default=100)

    global opts
    opts = parser.parse_args()

    in_file = open(opts.file)
    read_file_blast(in_file)
    
main()