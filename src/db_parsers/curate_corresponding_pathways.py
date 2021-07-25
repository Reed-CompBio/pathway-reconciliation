import sys

def main(INFILE,OUTFILE):

    out = open(OUTFILE,'w')
    header = None
    to_skip=[]
    with open(INFILE) as fin:
        for line in fin:
            if header == None:
                out.write(line)
                header = line.strip().split('\t')
                continue
            this_row = line.strip().split('\t')
            new_row = []
            pws_to_choose = display(this_row,header)
            this_name = input('NAME:')
            if len(this_name)=='':
                to_skip.append(this_row)
                continue
            new_row.append(this_name)
            for i in range(1,len(this_row)):
                if i in pws_to_choose:
                    pw_choice = this_row[i].split('|')
                    this_num = input('  CHOOSE FOR %s (%d-%d):' % (header[i],1,len(pw_choice)))
                    new_row.append(pw_choice[int(this_num)-1])
                else:
                    new_row.append(this_row[i])
            out.write('%s\n' % ('\t'.join(new_row)))
    out.close()
    print('skipped %d pathways:' % (len(to_skip)))
    for line in to_skip:
        print([n for n in line if n != 'NAME' and n != 'NaN'])
    print('DONE')
    return


def display(row,header):
    print('--')
    pws_to_choose = set()
    for i in range(1,len(row)):
        if row[i] == 'NaN':
            print('%s:' % (header[i]))
        elif '|' in row[i]:
            pws_to_choose.add(i)
            print('%s:' % (header[i]))
            n = 1
            for p in row[i].split('|'):
                print('  [%d] %s' % (n,p))
                n+=1
        else:
            print('%s: %s'  %(header[i],row[i]))
    print('--')
    return pws_to_choose

if __name__ == '__main__':
    #INFILE = '../../networks/dbs/corresponding-top-picks-ORIG.txt'
    #OUTFILE = '../../networks/dbs/corresponding-top-picks.txt'
    main(sys.argv[1],sys.argv[2])
