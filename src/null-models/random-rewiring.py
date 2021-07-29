#
# For every network in a given directory,
# this script generates random network that preserve the degree of each node
# example run:
# python random-rewiring.py ../../graphlets/dbs/ 10 ../../networks/null-models/random-rewiring/
#

'''
Original method (timed - only for Hedgehog 10):
real	0m7.477s
user	0m7.949s
sys	0m0.389s
'''

'''
New method (timed - only for Hedgehog 10):
real	0m2.887s
user	0m3.199s
sys	0m0.407s

New method (timed - only for Hedgehog 100):
real	0m18.456s
user	0m18.337s
sys	0m0.767s
'''

import os
import sys
import glob
import pandas as pd
import networkx as nx
import random

###################################
# create null model networks using configuration model or rewiring

def create_random_networks(ns: int, randdir: str, f: str, inputG: nx.Graph, method: str):

    #ns: number of randomized networks to be created (int)
    #randdir: output directory where the output is to be stored (string)
    #f: name of the original pathway file ending with '-edges.txt' (string)
    #inputG: Original Network (nx Graph)
    #
    #writes ns random network edgelists to the directory specified
    if method==1:

        k = [x for n,x in inputG.degree()]
        #print(k)
        for i in range(0,ns):

            Gr = nx.configuration_model(k)
            Gr = nx.Graph(Gr)
            Gr.remove_edges_from(nx.selfloop_edges(Gr))
            fullname = randdir+f+'-'+str(i)+'randomnetwork.txt'
            nx.write_edgelist(Gr,fullname,data=False,delimiter=" ")
            print(f+'-'+str(i)+'randomnetwork.txt')

    elif method==0:
        for i in range(0,ns):
            print('***',i)
            Gr=inputG
            nt=5*nx.number_of_edges(Gr)
            max_attempt = 50*nx.number_of_edges(Gr)

            x=0
            y=0
            while ((x < nt) and (y < max_attempt)):

                edges = list(Gr.edges)
                edge1 = random.choice(edges)
                edge2 = random.choice(edges)
                y=y+1
                if edge1!=edge2 and edge1[0]!=edge2[1] and edge2[0]!=edge1[1] and not Gr.has_edge(edge1[0], edge2[1]) and not Gr.has_edge(edge2[0], edge1[1]):

                    Gr.remove_edge(edge1[0], edge1[1])
                    Gr.remove_edge(edge2[0], edge2[1])

                    Gr.add_edge(edge1[0], edge2[1])
                    Gr.add_edge(edge2[0], edge1[1])
                    x=x+1
            fullname = randdir+f+'-'+str(i)+'randomnetwork.txt'
            if not os.path.isdir(os.path.dirname(fullname)):
                os.mkdir(os.path.dirname(fullname))
                print('Making Directory',os.path.dirname(fullname))

            nx.write_edgelist(Gr,fullname,data=False,delimiter=" ")
            print(x/nt, y/max_attempt, f+'-'+str(i)+'randomnetwork.txt')

def create_random_networks_2(ns: int, randdir: str, f: str, inputG: dict, num_edges: int, method: str):

    #ns: number of randomized networks to be created (int)
    #randdir: output directory where the output is to be stored (string)
    #f: name of the original pathway file ending with '-edges.txt' (string)
    #inputG: Original Network (nx Graph)
    #
    #writes ns random network edgelists to the directory specified
    if method==1:
        sys.exit('Method 1 for v2 is not implemented yet!')

    elif method==0:
        for i in range(0,ns):
            print('***',i)
            Gr = inputG.copy()
            nt=5*num_edges
            max_attempt = 50*num_edges

            x=0
            y=0
            while ((x < nt) and (y < max_attempt)):
                y=y+1

                # step 1: choose two nodes randomly.
                u,v = random.sample(Gr.keys(),2)

                # step 2: choose neighbors of those nodes that aren't each other.
                try:
                    u_neighbor = random.sample(Gr[u].difference(set([v])))
                    v_neighbor = random.sample(Gr[v].difference(set([u,u_neighbor])))
                except:
                    continue

                ## we now have two edges taht are guaranteed to be swappable. Swap!

                ## Remove (u,u_neighbor) and (v,v_neighbor)
                Gr[u].remove(u_neighbor)
                Gr[u_neighbor].remove(u)
                Gr[v].remove(v_neighbor)
                Gr[v_neighbor].remove(v)

                # Add (u,v_neighbor) and (v,u_neighbor)
                Gr[u].add(v_neighbor)
                Gr[v_neighbor].add(u)
                Gr[v].add(u_neighbor)
                Gr[u_neighbor].add(v)
                x=x+1

            fullname = randdir+f+'-'+str(i)+'randomnetwork.txt'
            if not os.path.isdir(os.path.dirname(fullname)):
                os.mkdir(os.path.dirname(fullname))
                print('Making Directory',os.path.dirname(fullname))

            out = open(fullname,'w')
            seen_edges = set()
            for u in Gr:
                for v in Gr[u]:
                    if (v,u) not in seen_edges:
                        out.write('%s %s\n' % (u,v))
                        seen_edges.add((u,v))
                        seen_edges.add((v,u))
            print(x/nt, y/max_attempt, f+'-'+str(i)+'randomnetwork.txt')

def main(argv):


    input_dir = argv[1] #r'../../graphlets/dbs/'
    iterations = int(argv[2])
    outdir = argv[3] #r'../../networks/null-models/' #randomnetwork.txt

    path_cp = r'../../networks/dbs/'
    correspondence_file='corresponding-top-picks-7pathways-6overlap-withcollapsed.txt'
    correspondence_file=path_cp+correspondence_file

    cp = pd.read_csv(correspondence_file,sep='\t',index_col=0)

    dbs = list(cp.columns.values)
    pathways = list(cp.index.values)

    for i in pathways:
        for j in dbs:
            if (pd.notnull(cp.loc[i,j])):
                f = input_dir+cp.loc[i,j][:-4]+'-network.txt'
                df=pd.read_csv(f, comment='#', header=None, sep='\s+')

                # new - make adjacency list. Note that we add (u,v) AND (v,u) in the list.
                G = {}
                num_edges = 0
                for index,row in df.iterrows():
                    if row[0] not in G:
                        G[row[0]] = set()
                    if row[1] not in G:
                        G[row[1]] = set()
                    G[row[0]].add(row[1])
                    G[row[1]].add(row[0])
                    num_edges+=1
                f = cp.loc[i,j][:-4]
                create_random_networks_2(iterations, outdir, f, G, num_edges, 0)

                # orig
                #G=nx.from_pandas_edgelist(df,0,1,create_using=nx.DiGraph()) #assuming the intial edgelist is directed
                #G=G.to_undirected()
                #G.remove_edges_from(nx.selfloop_edges(G))
                #f = cp.loc[i,j][:-4]
                #create_random_networks_2(iterations, outdir, f, G, 1)

        sys.exit() # for timing



if __name__ == "__main__":
    main(sys.argv)
