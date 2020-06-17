import argparse

def read_sam(sam):
    '''Read sam file'''

    for line in sam:
        
        if line.startswith('@'):
            print(line.rstrip())
        else:
            break


def main():
    parser = argparse.ArgumentParser('Extracts header from the sam file')
    parser.add_argument('-sam', '--file', help="sam file")
    
    global opts
    opts = parser.parse_args()
    
    sam = open(opts.file)
    read_sam(sam)
    
main()