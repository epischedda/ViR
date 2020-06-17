import argparse

class Hit():
    RG=''
    TEid=''
    score=''


def read_blast(blast):
    '''Read blast file'''

    HitDict={}
    c=0

    for line in blast:
        
        line=line.rstrip()
        row=line.split('\t')
        hit_ist=Hit()
        hit_ist.RG=row[0]
        hit_ist.TEid=row[1]
        hit_ist.aln_len=row[3]
        hit_ist.score=row[-1]
        temp_d=Hit()

        if int(hit_ist.aln_len) > int(opts.min_len):
            if hit_ist.RG not in HitDict.keys():
                HitDict[hit_ist.RG]=hit_ist
            else:
                temp_d=HitDict[hit_ist.RG]
                if float(temp_d.score) < float(hit_ist.score):
                    temp_d=hit_ist
                    HitDict[hit_ist.RG]=temp_d

    return HitDict


def update_groupInfo(groupInfo,HitDict):
    '''Read groupInfo file, add the information of the best match with TE and print a new output'''

    out_groups_info=open(opts.outputPath+'/Complete_Read_Groups_Info_final.txt','w')
    
    GroupsDict={}
    c=0
    
    for line in groupInfo:
        if line.startswith('READ_GROUP'):
            out_groups_info.write(line)
            continue
        elif line.startswith('Ungrouped'):
            out_groups_info.write(line)
            continue
        elif line.startswith('Shared'):
            out_groups_info.write(line)
            continue
        else:
            line=line.rstrip()
            row=line.split('\t')
            first_part='\t'.join(row[0:8])
            if row[0] in HitDict.keys():
                TE=HitDict[row[0]].TEid
            else:
                TE=''
            second_part='\t'.join(row[9:len(row)])
            
            out_groups_info.write(first_part+'\t'+str(TE)+'\t'+second_part+'\n')

    return HitDict

    
def main():
    parser = argparse.ArgumentParser('For each final read groups, add the information of the best match with TE in the Complete_Groups_Info.txt output file')
    parser.add_argument('-blast_file', '--blastfile', help="blast output file")
    parser.add_argument('-groupInfo', '--groupInfo', help="groupInfo file")
    parser.add_argument('-min_len', '--min_len', help="minimum alignment length required", default=100)
    parser.add_argument('-o', '--outputPath', help="absolute path for output files")

    global opts
    opts = parser.parse_args()
    
    blast = open(opts.blastfile)
    HitDict=read_blast(blast)
    
    groups = open(opts.groupInfo)
    update_groupInfo(groups,HitDict)

main()