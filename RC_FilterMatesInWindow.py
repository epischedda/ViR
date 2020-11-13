import argparse

def read_mates_in_window(mates_file):
    '''Read the output of bedtools closest between host and viral reads'''
    
    out_file=str(opts.outputPath)+'/List_FinalReads.txt'
    out=open(out_file,'w')
    ToRemoveReads=[]

    for line in mates_file:

        line=line.rstrip()
        row=line.split('\t')
        
        host=row[3].split('_bh_')[0]
        viral=row[9].split('_bh_')[0]

        if str(host)==str(viral):
            ToRemoveReads.append(host)

    ToRemoveReads = set(ToRemoveReads)
    for read in list(ToRemoveReads):
        out.write(str(read)+'\n')
    
    out.close()
  

def main():
    parser = argparse.ArgumentParser('Select pairs in which host and viral reads map in the reference genome within a window having dimension selected by the user')
    parser.add_argument('-i_mates_in_window', '--file_mate_in_window', help="output of bedtools closest between host and viral reads")
    parser.add_argument('-o', '--outputPath', help="output directory")
    
    global opts
    opts = parser.parse_args()

    mates_file = open(opts.file_mate_in_window)
    ReadsDict=read_mates_in_window(mates_file)
    
main()