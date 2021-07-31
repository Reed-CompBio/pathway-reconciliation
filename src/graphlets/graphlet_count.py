###################################
# Undirected Orbit & Graphlet frequency
#
# python graphlet_count.py /path/to/input/ /path/to/output_dir/
# input directory contains (pathway-edgelist).txt files.
###################################

import pandas as pd
import glob
import networkx as nx
import os
import argparse
import sys
import random

def gcount(orbitcount: str):
###################################
## graphlet counts from orbit count
## Input: orbit count (.ocount) file
##        generated from ORCA code
###################################
    print('running gcount')
    #read the .ocount file
    df=pd.read_table(orbitcount,sep=" ",header=None)

    # the different graphlet sizes (up to 5)
    # orca does not count larger than 5-node graphlets
    # check the number of columns in the .ocount files

    if len(df.columns)==1: #orbits 0-1
        d=[2]
        df=(df.iloc[:,0:1].sum(axis=0)/d).astype(int)
        df.to_csv(orbitcount[:-7]+'.gcount',sep=' ',header=False, index=True)

    elif len(df.columns)==4: #orbits 0-3
        d=[2,2,1,3]
        df=(df.iloc[:,0:4].sum(axis=0)/d).astype(int)
        df=df[[0,1,3]]
        df=df.reset_index(drop=True)
        df.to_csv(orbitcount[:-7]+'.gcount',sep=' ',header=False, index=True)

    elif len(df.columns)==15: #orbits 0-14
        d=[2,2,1,3,2,2,3,1,4,1,2,1,2,2,4]
        df=(df.iloc[:,0:15].sum(axis=0)/d).astype(int)
        df=df[[0,1,3,4,6,8,9,12,14]]
        df=df.reset_index(drop=True)
        df.to_csv(orbitcount[:-7]+'.gcount',sep=' ',header=False, index=True)

    elif len(df.columns)==73: #orbits 0-72
        d=[2,2,1,3,2,2,3,1,4,1,2,1,2,2,4, #0-14
        2,2,1,1,2,1,1,4,1,2,1,2,1,1,2,1,2,2,1, #-33
        5,1,1,2,1,1,2,1,1,4,1,1,1,1,2, #48
        3,2,2,1,2,3,2,1,3,1,2,2,1,1,2,2,1,2,2, #-67
        4,1,2,3,5] #-72
        df=(df.iloc[:,0:73].sum(axis=0)/d).astype(int)
        df=df[[0,1,3,4,6,8,9,12,14,15,18,22,24,27,31,34,35,39,43,45,49,51,54,56,59,62,65,68,70,72]]
        df=df.reset_index(drop=True)
        df.to_csv(orbitcount[:-7]+'.gcount',sep=' ',header=False, index=True)
    else:
        print('incorrect number of coulmns')

    print('wrote to {}'.format(orbitcount[:-7]+'.gcount'))
####################################

def gcount_directory(path_to_dir: str):
# gcount() for all .ocount files in a specified directory
# creates the corresponding .gcount files in the same directory

    files = glob.glob(path_to_dir + '*.ocount')
    files = [os.path.join(path_to_dir,x) for x in os.listdir(path_to_dir) if '.ocount' in x]
    print(files)
    for f in files:
        print('calling gcount on {}'.format(f))
        try:
            gcount(f)
            #remove the .ocount files if not needed anymore
            #os.system('rm '+f)
        except:
            print("{} has something with it".format(f))



def main(argv):

    # Specify the directories here:
    #input directory is required, this will be the output directory if one is not supplied
    path = argv[1]
    outdir = argv[2]

    #you can supply an output directory as well if you want
    #if len(argv) > 2:
    #    outdir = argv[2]
    #    if not os.path.exists(outdir):
    #        os.mkdir(outdir)

    #Pathways
    filenames = [f for f in glob.glob(path + '*.txt') if sum(1 for l in open(f)) > 1 and not 'pathway-names.txt' in f]

    for f in filenames:
        print('processing file: {}'.format(f))
        try:
            try:
                df=pd.read_csv(f,comment='#',header=None,sep='\s+')
            except:
                df=pd.read_csv(f,comment='#',header=None,sep='\r|\t')

            df = df[[0,1]]

            df = df[df[0] != '-']
            df =  df[df[0] != '--']
            df = df[df[1] != '-']
            df =  df[df[1] != '--']
            # create nx Graph

            G = nx.from_pandas_edgelist(df,0,1,create_using=nx.DiGraph()) #assuming the intial edgelist is directed

            G_new = G.to_undirected() # undirected
            G_new.remove_edges_from(nx.selfloop_edges(G_new))

            f=f.replace(path,'')
            outfile = f+'.id'
            # Store the original node IDs
            fullname = os.path.join(outdir, outfile)

            G_new = nx.convert_node_labels_to_integers(G_new, first_label=0, ordering='default', label_attribute='ID')
            pd.DataFrame(list(G_new.nodes('ID'))).to_csv(fullname,sep='\t',header=False,index=False)

            # Store the graph edgelist with numberic node IDs ([0, N-1])
            outfile=f[:-4]+'-network.txt'
            fullname = os.path.join(outdir, outfile)
            nx.write_edgelist(G_new,fullname,data=False,delimiter=" ")

        except Exception as e:
            print("There's something suspicious about {}".format(f))
            print(e)
    #create .orca .ocount files for the original pathway networks
    print('formatting files for orca and running orca')
    os.system('sh files_for_orca.sh '+outdir)
    os.system('sh execute_orca.sh 5 '+outdir) # for 4-node graphlets

    gcount_directory(outdir)

if __name__ == "__main__":
    main(sys.argv)
