import os
import sys

def main(INFILE,OUTFILE):
    out = open(OUTFILE,'w')
    with open(INFILE) as fin:
        for line in fin:
            if line[0] == '#': #header
                old_header = line.strip().split('\t')
                new_header = [h for h in old_header]
                for n in old_header:
                    if '_expanded' in n:
                        new_header.append(n.replace('expanded','collapsed'))
                out.write('%s\n' % ('\t'.join(new_header)))
                continue
            old_row = line.strip().split('\t')
            new_row = [n for n in old_row]
            for n in old_row:
                if '_expanded' in n:
                    newfile = n.replace('expanded','collapsed')
                    if os.path.isfile('../../networks/dbs/'+newfile):
                        new_row.append(newfile)
                    else:
                        print('WARNING: %s does not exist.' % (newfile))
                        new_row.append('NaN')
            out.write('%s\n' % ('\t'.join(new_row)))
    out.close()
    print('done writing ',OUTFILE)
    return


if __name__ == '__main__':
    #INFILE = '../../networks/dbs/corresponding-top-picks.txt'
    #OUTFILE = '../../networks/dbs/corresponding-top-picks-withcollapsed.txt'
    main(sys.argv[1],sys.argv[2])
