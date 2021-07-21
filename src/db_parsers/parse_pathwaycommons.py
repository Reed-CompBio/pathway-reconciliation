import sys
import os
import argparse
import itertools
import glob

def main(args):
    infile = args.infile
    db_name = infile.split('.')[1]
    print('INFILE:',infile)
    print('DB NAME:',db_name)

    outdir = args.outdir+db_name
    if not os.path.isdir(outdir):
        print('making output directory %s...' % (outdir))
        os.makedirs(outdir)

    if db_name == 'All':
        outfile = outdir+'all_pathway_commons.txt'
        proteins = get_uniprot(infile,outdir+'/uniprot.txt')
        pathways = read_interactions(infile,proteins,all_file=True)
        write_file(pathways['all'],outfile)
    else:
        proteins = get_uniprot(infile,None)
        pathways = read_interactions(infile,proteins)
        out = open(outdir+'/pathway-names.txt','w')
        for p in pathways:
            if len(pathways[p]) < args.thres:
                #print('IGNORING %s: %d interactions.' % (p,len(pathways[p])))
                continue
            pwname = p.replace(' ','_').replace('/','-').replace('(','~').replace(')','~')+'.txt'
            out.write('%s\t%d\t%s\n' % (pwname,len(pathways[p]),p))
            write_file(pathways[p],outdir+'/'+pwname)
        out.close()

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
