import glob
import os
import csv

pp_relationships = ['controls-state-change-of','controls-phosphorylation-of',
    'controls-transport-of','controls-expression-of','catalysis-precedes',
    'in-complex-with','interacts-with']
def main():

    outdir = '../../networks/dbs/pathbank/'
    if not os.path.isdir(outdir):
        print('making output directory %s...' % (outdir))
        os.makedirs(outdir)

    mapper = get_uniprot_names()

    #pathway_names = get_pathway_names()

    files = glob.glob('infiles/pathbank/*.extended-sif')
    pathways = read_interactions(files,mapper)

    out = open(outdir+'/pathway-names.txt','w')
    num_small = 0
    tot=0
    #written_edgesets = set()
    for p in pathways:
        if len(pathways[p]) < 10:
            print('WARNING! pathway %s has %d interactions. Ignoring.' % (p,len(pathways[p])))
            num_small +=1
            continue
        pwname = p.replace(' ','_').replace('/','-').replace('(','~').replace(')','~')+'.txt'
        out.write('%s\t%d\t%s\n' % (pwname,len(pathways[p]),p))
        write_file(pathways[p],outdir+'/'+pwname)
        tot+=1
    out.close()
    print('%d too small. %d TOTAL parsed.' % (num_small,tot))
    return

def get_pathway_names():
    all_pathway_names = {}
    with open('infiles/pathbank/pathbank_pathways.csv') as fin:
        reader = csv.reader(fin, delimiter=',', quotechar='"')
        for row in reader:
            if 'SMPDB' in row[0]:
                continue # skip header
            if len(row) > 2:
                print(row[1],row[2])
                all_pathway_names[row[1]] = row[2]
    return all_pathway_names

'''
Parsed undirected edges from SIF.
Sorts each edge to ensure that only one of (u,v) or (v,u) are added.
Only considers edges where both nodes are in the proteins dictionary.
Does not handle the MEDIATOR_IDS.
'''
def read_interactions(files,mapper):
    pathways = {} # pathway to edge tuple
    print(len(files))
    for infile in files:
        num_missed = 0
        num_processed = 0
        p = None
        with open(infile) as fin:
            print(infile)
            for line in fin:
                if 'PARTICIPANT_TYPE' in line:
                    # done reading interactions; return.
                    break
                if 'INTERACTION_TYPE' in line:
                    # skip header
                    continue

                row = line.strip().split('\t')
                if len(row)>5:
                    p = row[5]
                    if p not in pathways:
                        pathways[p] = set()

                if len(row)> 2 and not any([row[1]==rel for rel in pp_relationships]):
                    # must be a protein-protein relationship.
                    continue

                if row[0] not in mapper or row[2] not in mapper:
                    num_missed +=1
                    continue

                num_processed+=1
                edge = tuple(sorted([mapper[row[0]],mapper[row[2]]]))
                pathways[p].add(edge)
            if p != None:
                print('  %s --> %d interactions; %d processed & %d skipped' % (p,len(pathways[p]),num_processed,num_missed))

    ## this should never happen (we should return when we hit PARTICIPANT_TYPE).
    ## return None in this case, it's an error.
    return pathways

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
