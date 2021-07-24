import glob
import os

pp_relationships = ['controls-state-change-of','controls-phosphorylation-of',
    'controls-transport-of','controls-expression-of','catalysis-precedes',
    'in-complex-with','interacts-with']
def main():

    outdir = '../../networks/dbs/reactome/'
    if not os.path.isdir(outdir):
        print('making output directory %s...' % (outdir))
        os.makedirs(outdir)

    mapper = get_uniprot_names()
    #pathways = read_interactions('infiles/reactome-wnt.extended-sif',mapper)
    pathways = read_interactions('infiles/reactome.extended-sif',mapper)
    files = glob.glob('infiles/netpath/*')
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

'''
Parsed undirected edges from SIF.
Sorts each edge to ensure that only one of (u,v) or (v,u) are added.
Only considers edges where both nodes are in the proteins dictionary.
Does not handle the MEDIATOR_IDS.
'''
def read_interactions(infile,mapper,all_file=False):
    pathways = {} # pathway to edge tuple
    if all_file:
        pathways['all'] = set()
    DELIM = ';'
    num_missed = 0
    with open(infile) as fin:
        for line in fin:
            if 'PARTICIPANT_TYPE' in line:
                # done reading interactions; return.
                return pathways
            if 'INTERACTION_TYPE' in line:
                # skip header
                continue

            row = line.strip().split('\t')
            if len(row)> 2 and not any([row[1]==rel for rel in pp_relationships]):
                # must be a protein-protein relationship.
                continue
            if row[0] not in mapper or row[2] not in mapper:
                num_missed +=1
                continue

            edge = tuple(sorted([mapper[row[0]],mapper[row[2]]]))

            if all_file:
                pathways['all'].add(edge)
                continue

            for p in row[5].split(';'):
                p = p.strip()
                if p == '':
                    continue
                if p not in pathways:
                    pathways[p] = set()
                pathways[p].add(edge)

    ## this should never happen (we should return when we hit PARTICIPANT_TYPE).
    ## return None in this case, it's an error.
    return None

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
