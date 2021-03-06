import sys
import os
import argparse
import itertools
import glob

# words to ignore for INOH pathways (these are redundant)
words_to_ignore = ['C. elegans','Drosophlia','Mammalian','Xenopus']

def main(args):
    infile = args.infile
    db_name = infile.split('.')[1]
    print('INFILE:',infile)
    print('DB NAME:',db_name)



    if db_name == 'All':
        outdir = args.outdir
        if not os.path.isdir(outdir):
            print('making output directory %s...' % (outdir))
            os.makedirs(outdir)

        outfile = outdir+'All_Pathway_Commons.txt'
        proteins = get_uniprot(infile,outdir+'Pathway_Commons_Uniprot_IDs.txt')
        pathways = read_interactions(infile,proteins,all_file=True)
        write_file(pathways['all'],outfile)
    else:
        outdir = args.outdir+db_name
        if not os.path.isdir(outdir):
            print('making output directory %s...' % (outdir))
            os.makedirs(outdir)

        num_small = 0
        num_nonhuman= 0
        num_identical=0
        tot=0
        proteins = get_uniprot(infile,None)
        pathways = read_interactions(infile,proteins)
        out = open(outdir+'/pathway-names.txt','w')
        written_edgesets = set()
        for p in pathways:
            if len(pathways[p]) < args.thres:
                print('WARNING! pathway %s has %d interactions. Ignoring.' % (p,len(pathways[p])))
                num_small +=1
                continue
            if db_name=='inoh' and any([w in p for w in words_to_ignore]):
                print('WARNING! Non-human pathway %s. Ignoring.' % (p))
                num_nonhuman+=1
                continue
            if frozenset(pathways[p]) in written_edgesets:
                print('WARNING! Identical pathway %s. Ignoring.' % (p))
                num_identical+=1
                continue
            pwname = p.replace(' ','_').replace('/','-').replace('(','~').replace(')','~')+'.txt'
            out.write('%s\t%d\t%s\n' % (pwname,len(pathways[p]),p))
            write_file(pathways[p],outdir+'/'+pwname)
            tot+=1
            written_edgesets.add(frozenset(pathways[p]))
        out.close()
        print('%d too small, %d nonhuman, %d identical. %d TOTAL parsed.' % (num_small,num_nonhuman,num_identical,tot))

    return


'''
Get uniprot-mapped protein names. Only use these in pathways.
Input: "All" file
Output: Dictionary of (name,uniprotID) pairs.
'''
def get_uniprot(infile,outfile):
    proteins = {}
    with open(infile) as fin:
        for line in fin:
            if 'uniprot knowledgebase:' in line:
                row = line.strip().split('\t')
                proteins[row[0]] = row[3].split(':')[-1]
    print('%d uniprot' % (len(proteins)))
    if outfile:
        fout = open(outfile,'w')
        for n in proteins:
            fout.write('%s\t%s\n' % (n,proteins[n]))
        fout.close()
    return proteins

'''
Parsed undirected edges from SIF.
Sorts each edge to ensure that only one of (u,v) or (v,u) are added.
Only considers edges where both nodes are in the proteins dictionary.
Does not handle the MEDIATOR_IDS.
'''
def read_interactions(infile,proteins,all_file=False):
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
            if row[0] not in proteins or row[2] not in proteins:
                num_missed +=1
                continue

            edge = tuple(sorted([row[0],row[2]]))

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

def write_file(edges,outfile):
    out = open(outfile,'w')
    for e in edges:
        out.write('%s\t%s\n' % (e[0],e[1]))
    out.close()
    print('  wrote %d lines to %s' % (len(edges),outfile))
    return

def parse_arguments():
    """
    Argument Parser for parse_pc.py.

    Returns
    -----------
    ArgumentParser object

    """
    parser = argparse.ArgumentParser('PathwayCommons Parser.  Converts pathways into undirected graphs with UniProtKB identifiers.')
    parser.add_argument('-i','--infile',help='input file (PathwayCommons11.*.hgnc.txt). Required')
    parser.add_argument('-o','--outdir',
        help='outfile directory. Default=out.',default='out/')
    parser.add_argument('-t','--thres',type=int,default=10,
        help='Do not write pathways with fewer than THRES edges. Default 10.')
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        print('making output directory %s...' % (args.outdir))
        os.makedirs(args.outdir)
    return args

if __name__ == '__main__':
    main(parse_arguments())
