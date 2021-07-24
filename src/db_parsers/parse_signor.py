import sys
import os
import argparse
import itertools
import glob
import json
import ndex2
import itertools

def main():

    mapper = get_uniprot_names()

    outdir_collapsed = '../../networks/dbs/signor_collapsed'
    if not os.path.isdir(outdir_collapsed):
        print('making output directory %s...' % (outdir_collapsed))
        os.makedirs(outdir_collapsed)

    outdir_expanded = '../../networks/dbs/signor_expanded'
    if not os.path.isdir(outdir_expanded):
        print('making output directory %s...' % (outdir_expanded))
        os.makedirs(outdir_expanded)


    families,family_ids = get_families(mapper)
    complexes = get_complexes(mapper,family_ids)
    pathways_collapsed,pathways_expanded = read_interactions('infiles/signor/*.cx',families,complexes)

    num_small = 0
    tot=0
    out = open('../../networks/dbs/signor_collapsed/pathway-names.txt','w')
    for p in pathways_collapsed:
        if len(pathways_collapsed[p]) < 10:
            print('WARNING! pathway %s has %d interactions. Ignoring.' % (p,len(pathways_collapsed[p])))
            num_small +=1
            continue
        pwname = p.replace(' ','_').replace('/','-').replace('(','~').replace(')','~')+'.txt'
        out.write('%s\t%d\t%s\n' % (pwname,len(pathways_collapsed[p]),p))
        write_file(pathways_collapsed[p],'../../networks/dbs/signor_collapsed/'+pwname)
        tot+=1
    out.close()
    print('signor_collapsed: %d too small,. %d TOTAL parsed.' % (num_small,tot))


    tot = 0
    num_small = 0
    out = open('../../networks/dbs/signor_expanded/pathway-names.txt','w')
    for p in pathways_expanded:
        if len(pathways_expanded[p]) < 10:
            print('WARNING! pathway %s has %d interactions. Ignoring.' % (p,len(pathways_expanded[p])))
            num_small +=1
            continue
        pwname = p.replace(' ','_').replace('/','-').replace('(','~').replace(')','~')+'.txt'
        out.write('%s\t%d\t%s\n' % (pwname,len(pathways_expanded[p]),p))
        write_file(pathways_expanded[p],'../../networks/dbs/signor_expanded/'+pwname)
        tot+=1
    out.close()
    print('signor_expanded: %d too small,. %d TOTAL parsed.' % (num_small,tot))

    return

def get_families(mapping):
    families = {}
    family_ids = {} # for mapping complexes
    with open('infiles/signor/SIGNOR_PF.csv') as fin:
        for line in fin:
            if 'SIGNOR ID' in line:
                continue # skip header
            row = line.strip().split(';')
            family_name = row[1]
            protein_ids = row[2].replace('"','').replace(' ','').split(',')
            protein_names = [mapping[p] for p in protein_ids]
            families[family_name] = set(protein_names)
            family_ids[row[0]] = set(protein_names)
    print('%d family names' % (len(families)))
    return families,family_ids


def get_complexes(mapping,families):
    complexes = {}
    complex_ids = {}
    with open('infiles/signor/SIGNOR_complexes.csv') as fin:
        for line in fin:
            if 'SIGNOR ID' in line:
                continue # skip header
            row = line.strip().split(';')
            complex_name = row[1]
            protein_ids = row[2].replace('"','').replace(' ','').split(',')
            protein_names = set()
            for p in protein_ids:
                if p in mapping:
                    protein_names.add(mapping[p])
                elif p in families:
                    protein_names.update(families[p])
                elif 'SIGNOR-C' in p: # another complex
                    protein_names.add(p)
                # ignore all other itms (mismapped uniprot IDs, empty spaces, chem names).
            complexes[complex_name] = set(protein_names)
            complex_ids[row[0]]=set(protein_names)


    # go back over complexes and replace 'SIGNOR-C' complexes with members.
    for c in complexes:
        to_remove = set()
        to_add = set()
        for member in complexes[c]:
            if 'SIGNOR-C' in member:
                to_add.update(complex_ids[member])
                to_remove.add(member)
        #if len(to_remove)>0:
    #        print('replacing',to_remove,'with',to_add)
        complexes[c].update(to_add) ## add all members
        for n in to_remove:
            complexes[c].remove(n) # removing complexID
    print('%d complexes names' % (len(complexes)))
    return complexes

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


'''
Parsed undirected edges, both collapsed and expanded.
Sorts each edge to ensure that only one of (u,v) or (v,u) are added.
'''
def read_interactions(pattern,families,complexes):
    pathways = {} # pathway to edge tuple
    c_edges = {}
    e_edges = {}

    for f in glob.glob(pattern):
        print(f)
        # Convert downloaded network to NiceCXNetwork object
        net_cx = ndex2.create_nice_cx_from_file(f)
        name = net_cx.get_name()
        # Display information about network and output 1st 100 characters of CX
        print('Name: ', net_cx.get_name())
        print('Number of nodes: ', str(len(list(net_cx.get_nodes()))))
        #print(net_cx.get_nodes())
        print('Number of edges: ' ,str(len(list(net_cx.get_edges()))))
        #print(net_cx.get_edges())

        c_edges[name] = set()
        e_edges[name] = set()

        # Create SIF network
        #nodes
        #0 {'@id': 0, 'n': 'DVL1', 'r': 'uniprot:O14640'}
        nodes = {}

        for k,v in net_cx.get_nodes():
            if 'SIGNOR-ST' in v['r'] or 'SIGNOR-PH' in v['r']:
                print('ignoring %s (%s)' % (v['n'],v['r']))
                continue # ignore phenotype and stimuli information
            nodes[k]={'name':v['n'],'type':v['r']}
        print('%d nodes' % (len(nodes)))

        parsed_complexes = set()
        for k,v in net_cx.get_edges():
            if v['s'] not in nodes or v['t'] not in nodes:
                continue
            node1 = nodes[v['s']]['name']
            node2 = nodes[v['t']]['name']
            edge = tuple(sorted([node1,node2]))
            c_edges[name].add(edge)

            if node1 in families:
                node1_type = 'FAMILY'
            elif node1 in complexes:
                node1_type = 'COMPLEX'
            else:
                node1_type = 'PROTEIN'

            if node2 in families:
                node2_type = 'FAMILY'
            elif node2 in complexes:
                node2_type = 'COMPLEX'
            else:
                node2_type = 'PROTEIN'

            if node1_type == 'PROTEIN' and node2_type == 'PROTEIN':
                # collapsed edge is same as expanded edge.
                e_edges[name].add(edge)
                continue

            node1_to_connect = set()
            if node1_type == 'PROTEIN':
                node1_to_connect.add(node1)
            elif node1_type == 'COMPLEX':
                node1_to_connect.update(complexes[node1])
                if node1 not in parsed_complexes: # add all v. all interactions
                    for n1,n2 in itertools.combinations(node1_to_connect,2):
                        edge = tuple(sorted([n1,n2]))
                        e_edges[name].add(edge)
                    parsed_complexes.add(node1)
            elif node1_type == 'FAMILY':
                node1_to_connect.update(families[node1])

            node2_to_connect = set()
            if node2_type == 'PROTEIN':
                node2_to_connect.add(node2)
            elif node2_type == 'COMPLEX':
                node2_to_connect.update(complexes[node2])
                if node2 not in parsed_complexes: # add all v. all interactions
                    for n1,n2 in itertools.combinations(node2_to_connect,2):
                        edge = tuple(sorted([n1,n2]))
                        e_edges[name].add(edge)
                    parsed_complexes.add(node2)
            elif node2_type == 'FAMILY':
                node2_to_connect.update(families[node2])
            #
            for n1 in node1_to_connect:
                for n2 in node2_to_connect:
                    edge = tuple(sorted([n1,n2]))
                    e_edges[name].add(edge)
        print('%d collapsed edges' % (len(c_edges[name])))
        print('%d expanded edges' % (len(e_edges[name])))
        print('--'*10)


    ## this should never happen (we should return when we hit PARTICIPANT_TYPE).
    ## return None in this case, it's an error.
    return c_edges,e_edges

def write_file(edges,outfile):
    out = open(outfile,'w')
    for e in edges:
        out.write('%s\t%s\n' % (e[0],e[1]))
    out.close()
    print('  wrote %d lines to %s' % (len(edges),outfile))
    return

if __name__ == '__main__':
    main()
