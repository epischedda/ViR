import argparse

class HIT:
    scaffold=''
    start=0
    end=0
    count=0
    reads=[]
    seq_scaffold=''
    seq_start=0
    seq_end=0
    seq_name=''
    seq_distance=0
    pass


def read_file_merge(txt):
    '''Read output bedtools merge'''
    out_file_tot=open(opts.outputFile,'w')
    out_file_tot.write('scaffold\tstart\tend\tid\tcount\treads\t'+str(opts.seq_type)+'_name'+'\t'+str(opts.seq_type)+'_start'+'\t'+str(opts.seq_type)+'_end'+'\t'+str(opts.seq_type)+'_distance'+'\n')
    
    hit_db={}

    for line in txt:
        
        reads_list=[]
        line=line.rstrip()
        parts=line.split('\t')
        
        hit_istance=HIT()
        hit_istance.scaffold=parts[0]
        hit_istance.start=parts[1]
        hit_istance.end=parts[2]
        hit_istance.count=len(parts[3].split(','))
        r_list=parts[3].split(',')
        hit_istance.seq_scaffold=parts[4]
        hit_istance.seq_start=parts[5]
        hit_istance.seq_end=parts[6]
        hit_istance.seq_name=parts[7]
        hit_istance.seq_distance=parts[8]
        
        for read in r_list:
            name=read.split('_')[0]
            reads_list.append(name)
        
        real_list = set(sorted(reads_list))
        hit_istance.reads=real_list
        hit_id=str(hit_istance.scaffold)+':'+str(hit_istance.start)+'-'+str(hit_istance.end)
        
        hit_db[hit_id]=hit_istance
        
        new_line=str(hit_istance.scaffold)+'\t'+str(hit_istance.start)+'\t'+str(hit_istance.end)+'\t'+str(hit_id)+'\t'+str(len(real_list))+'\t'+','.join(real_list)
        
        if int(hit_istance.seq_distance) <= int(opts.max_dist) and int(hit_istance.seq_distance) != -1:
            nirv_info='\t'+str(hit_istance.seq_name)+'\t'+str(hit_istance.seq_start)+'\t'+str(hit_istance.seq_end)+'\t'+str(hit_istance.seq_distance)
            
        else:
            nirv_info='\t.\t.\t.\t.'
            
        out_file_tot.write(new_line+nirv_info+'\n')
            
    out_file_tot.close()

    return hit_db
  
    
def main():
    parser = argparse.ArgumentParser('Filter bedtools closest output by distance as selected by the user')
    parser.add_argument('-i', '--file', help="bedtools merge output input file")
    parser.add_argument('-o', '--outputFile', help="output path file")
    parser.add_argument('-d', '--max_dist', help="max distance between the coordinate of the equivalent region and an EVE annotated", default=10000)
    parser.add_argument('-seq_type', '--seq_type', help="type of sequences represented near each blast result shown")
    
    global opts
    opts = parser.parse_args()

    in_file = open(opts.file)
    read_file_merge(in_file)
    
main()