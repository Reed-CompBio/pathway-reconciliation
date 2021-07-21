import glob
import os

def main():
    pathway_names = get_pathway_names('infiles/kegg/_pathway_list.txt')
    mapper = get_uniprot_names()
    process('collapsed',pathway_names,None)
    process('expanded',pathway_names,mapper)

def process(ptype,pathway_names,mapper): # collapsed or expanded
    # assumes files are in infiles/netpath.
    outdir = '../../networks/dbs/kegg_%s/' % (ptype)
    if not os.path.isdir(outdir):
        print('making output directory %s...' % (outdir))
        os.makedirs(outdir)

    files = glob.glob('infiles/kegg/*%s-edges.txt'% (ptype))
    out = open(outdir+'/pathway-names.txt','w')
    for f in files:
        p = f.replace('infiles/kegg/','').replace('-%s-edges.txt' % (ptype),'')
        pname = pathway_names[p]
        #print(pname)
        pwname = pname.replace(' ','_').replace('/','-').replace('(','~').replace(')','~')+'.txt'
        edges = set()
        missed = 0
        with open(f) as fin:
            for line in fin:
                if line[0] == '#':
                    continue
                row = line.strip().split('\t')
                if mapper and row[0] in mapper and row[1] in mapper:
                    edge = tuple(sorted([mapper[row[0]],mapper[row[1]]]))
                    edges.add(edge)
                elif mapper == None and row[0] != '' and row[1] != '':
                    edge = tuple(sorted([row[0],row[1]]))
                    edges.add(edge)
                else:
                    missed+=1
        out.write('%s\t%d\t%s\n'% (pwname,len(edges),pname))
        outfile = outdir+'/'+pwname
        if missed > 0:
            print('  %d missed' % (missed))
        write_file(edges,outfile)
    out.close()

def get_pathway_names(infile):
    pnames = {}
    with open(infile) as fin:
        for line in fin:
            row =line.strip().split('\t')
            if len(row) != 2:
                continue
            row[0] = row[0].replace('path:','')
            row[1] = row[1].replace(' - Homo sapiens (human)','')
            pnames[row[0]] = row[1]
    return pnames

def get_uniprot_names():
    mapper = {}
    with open('mapping.txt') as fin:
        for line in fin:
            row = line.strip().split('\t')
            if len(row) != 2:
                continue
            mapper[row[1]]=row[0]
    print('%d uniprot-to-common names' % (len(mapper)))
    return mapper

def write_file(edges,outfile):
    out = open(outfile,'w')
    for e in edges:
        out.write('%s\t%s\n' % (e[0],e[1]))
    out.close()
    print('  wrote %d lines to %s' % (len(edges),outfile))
    return

if __name__ == '__main__':
    main()
