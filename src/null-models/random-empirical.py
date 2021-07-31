
#
# This script will generate random empirical networks to match a database
# Call it by running python3 random-empirical.py /path/to/interactome /path/to/db number-of-iterations /path/to/write
# 


import os
import sys
import pandas as pd
import networkx as nx
from math import log as math_log
import random
import re

## copied from PerfectLinker
def df_to_graph(fname: str,verbose=True,weighted=False):
        """
        :fname   path/to/dataframe
        :returns nx graph
        """
        if weighted:
                df = pd.read_csv(fname,sep='\t').take([0,1,2],axis=1)
        else:
                df = pd.read_csv(fname,sep='\t').take([0,1],axis=1)
        #df = pd.read_csv(fname,sep='\t',skiprows=0,names=['head','tail','weight']).take([0,1,2],axis=1)
        #df = df.drop(0)
        ff = df
        edges = [tuple(x) for x in df.values]
        #print(edges)
        vertices = list(df.stack())
        GR = nx.DiGraph()
        for v in vertices:
                GR.add_node(v)
        for e in edges:
                if weighted:
                        u,v,w = e
                        GR.add_edge(u,v,weight=float(w),cost=-math_log(max([0.000000001, w]))/math_log(10))
                else:
                        u,v = e
                        GR.add_edge(u,v)
        if verbose:
                print('length of edge set: {}'.format(len(set(edges))))
                print('number of edges in GR: {}'.format(len(list(GR.edges))))
        GR.remove_edges_from(nx.selfloop_edges(GR))
        return GR



## write the network
def write_output(edges,outfile,verbose=False):
        out = open(outfile,'w')
        out.write('#tail\thead\n')
        for u,v in edges:
                out.write('%s\t%s\n' % (u,v,))
        out.close()
        if verbose:
                print('wrote to %s' % (outfile))
        return


def construct_graph(G,s):
    H = set()
    at = random.choice(list(G.nodes))
    while len(H) < s:
        poss = list(nx.all_neighbors(G,at))
        nxt  = random.choice(poss)
        H.add(nxt)
        at = nxt
    return G.subgraph(H)


def construct_fancy_graph(G,s):
    V,E = s
    H = set()
    at = random.choice(list(G.nodes))
    while len(H) < V:
        poss = list(nx.all_neighbors(G,at))
        nxt  = random.choice(poss)
        H.add(nxt)
        at = nxt
    #okay, now we have the induced subgraph
    sub = G.subgraph(H).copy()
    #next, let's remove edges randomly
    itrs = 0
    while len(list(sub.edges)) > E:
        if itrs == 100:
            return sub
        at = random.choice(list(sub.edges))
        u,v = at
        if sub.degree[u] > 1 and sub.degree[v] > 1:
                sub.remove_edge(u,v)
                itrs = 0
        else:
            #print(itrs)
            itrs += 1

    return sub


def construct_graph_edges(G,s):
    H = []
    at = random.choice(list(G.nodes))
    while len(H) < s:
        poss = list(nx.all_neighbors(G,at))
        nxt  = random.choice(poss)
        H.append((at,nxt))
        at = nxt
    return G.edge_subgraph(H)

def construct_graph_restarts(G,s):
    H = set()
    p=0.0
    at = random.choice(list(G.nodes))
    strt = at 
    while len(H) < s:
        poss = list(nx.all_neighbors(G,at))
        restart = random.random() < p
        if restart:
            nxt = strt
        else:
            nxt  = random.choice(poss)
        H.add(nxt)
        at = nxt
    return G.subgraph(H)

def fetch_sizes(pathname):
    idfiles = [os.path.join(pathname,x) for x in os.listdir(pathname) if '.id' in x]
    d = dict()
    for i in idfiles:
        #nm = '-'.join(i.split('/')[-1].split('-')[:-1])
        #just let the name by the file
        nm = i.split('/')[-1].split('.')[0]
        with open(i,'r') as f:
            sz = len(f.read().splitlines())
        d[nm] = sz
    return d


def fancy_sizes(pathname):
    """
    :this version looks at the networks, rather than the id files. It controls for nodes and edges
    """
    netfiles = [os.path.join(pathname,x) for x in os.listdir(pathname) if re.search('-network.txt$',x) and not 'pathway-names' in x] 
    d = dict()
    for n in netfiles:
        #nm = '-'.join(i.split('/')[-1].split('-')[:-1])
        #just let the name by the file
        nm = n.split('/')[-1].split('.')[0]
        with open(n,'r') as f:
            E = len(f.read().splitlines())
        with open(n,'r') as f:
            V = len(set(f.read().split()))
        d[nm] = (V,E)
    return d

def gen_random_networks(G,S,it,out):
    for net in S:
        print(net)
        size = S[net]
        for i in range(it):
            #H = construct_graph(G,size)
            H = construct_fancy_graph(G,size)
            outpath = os.path.join(out,'{}-RandomWalkerInduced-{}-edges.txt'.format(net,i))
            write_output(H.edges,outpath)


def main(argv):
    #G = df_to_graph(argv[1])
    G = nx.read_edgelist(argv[1])
    #sizes = fetch_sizes(argv[2])
    sizes = fancy_sizes(argv[2])
    print(sizes)
    iterations = int(argv[3])
    print(iterations)
    out = argv[4]
    gen_random_networks(G,sizes,iterations,out)
if __name__ == "__main__":
    main(sys.argv)

