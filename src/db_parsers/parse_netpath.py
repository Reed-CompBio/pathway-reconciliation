import glob
import os

def main():
    # assumes files are in infiles/netpath.
    outdir = '../../networks/netpath/'
    if not os.path.isdir(outdir):
        print('making output directory %s...' % (outdir))
        os.makedirs(outdir)

    files = glob.glob('infiles/netpath/*')
    out = open(outdir+'/pathway-names.txt','w')
    for f in files:
        pname = f.replace('infiles/netpath/','').replace('-edges.txt','')
        #print(pname)
        out.write('%s-edges.txt\t%s\n'% (pname,pname))
        edges = set()
        with open(f) as fin:
            for line in fin:
                if line[0] == '#':
                    continue
                row = line.strip().split('\t')
                edge = tuple(sorted([row[6],row[7]]))
                edges.add(edge)
        outfile = outdir+pname+'.txt'
        write_file(edges,outfile)
    out.close()

def write_file(edges,outfile):
    out = open(outfile,'w')
    for e in edges:
        out.write('%s\t%s\n' % (e[0],e[1]))
    out.close()
    print('  wrote %d lines to %s' % (len(edges),outfile))
    return

if __name__ == '__main__':
    main()
