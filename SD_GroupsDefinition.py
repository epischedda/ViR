import argparse

class HIT:
    scaffold=''
    start=0
    end=0
    count=0
    reads=[]
    pass


def read_file_merge(txt):
    '''Read output of bedtools merge and print output files'''
    out_file_bed=open(opts.outputBed,'w')
    out_file_reginfo=open(opts.outputRegInfo,'w')
    out_file_reginfo.write('scaffold start end id count reads eve_name eve_start eve_end eve_distance piwi_name piwi_start piwi_end piwi_distance\n')
    hit_db={}

    for line in txt:
        
        reads_list=[]
        line=line.rstrip()
        parts=line.split('\t')
        
        hit_istance=HIT()
        hit_istance.scaffold=parts[0]
        hit_istance.start=parts[1]
        hit_istance.end=parts[2]
        hit_istance.count=parts[3]
        r_list=parts[4].split(',')
        
        for read in r_list:
            name=read.split('_')[0]
            reads_list.append(name)
        
        real_list = set(sorted(reads_list))
        hit_istance.reads=real_list
        hit_id=str(hit_istance.scaffold)+':'+str(hit_istance.start)+'-'+str(hit_istance.end)
        
        hit_db[hit_id]=hit_istance
        
        if int(len(real_list)) >= int(opts.minReads):
            line_bed=str(hit_istance.scaffold)+'\t'+str(hit_istance.start)+'\t'+str(hit_istance.end)+'\t'+','.join(real_list)+'\n'
            line_RegInfo=str(hit_istance.scaffold)+'\t'+str(hit_istance.start)+'\t'+str(hit_istance.end)+'\t'+str(hit_id)+'\t'+str(len(real_list))+'\t'+','.join(real_list)+'\t.\t.\t.\t.\t.\t.\t.\t.\n'
            out_file_bed.write(line_bed)
            out_file_reginfo.write(line_RegInfo)
        
    out_file_bed.close()
    out_file_reginfo.close()

    return hit_db

    
def main():
    parser = argparse.ArgumentParser('Filter out regions in which the number of reads ar less than the threshold')
    parser.add_argument('-i', '--file', help="bedtools merge output input file")
    parser.add_argument('-o_bed', '--outputBed', help="output path file of the bed file")
    parser.add_argument('-o_region_info', '--outputRegInfo', help="output path file of the region info file")
    parser.add_argument('-minReads', '--minReads', help="minimum number of reads in region", default=2)
    
    global opts
    opts = parser.parse_args()

    in_file = open(opts.file)
    read_file_merge(in_file)
    
main()