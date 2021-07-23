# python rho_coeff.py /path/to/orca_dir/ /path/to/output_dir/ /path/to/network/
# the third argument is optional - to be used only if -network.txt files are not in orca 
#
import pandas as pd 
import networkx as nx
import matplotlib.pyplot as plt
import glob
import sys
from scipy import stats

def compute_rho(path_to_network: str, path_to_orca: str, path_to_outdir: str):
    
    filenames = glob.glob(path_to_orca + '*.ocount')

    for f in filenames:

        df = pd.read_csv(f,sep=" ",header=None)

        f_network = f.split('/')[-1][:-12]
        print(f_network)
        fullpath = path_to_network+f_network
        df_network = pd.read_csv(fullpath,sep=" ",header=None)
        df_graph = nx.from_pandas_edgelist(df_network,0,1)

        P=(df>0).astype(int)

        if (P[0].sum()>0):
            rho1 = 1 - (2*(P[0].sum())/(df[0].sum()))
            #print(f, P[0].sum(),df[0].sum())
        else:
            rho1=0.0

        if ((P[1]*(1-P[3])).sum() > 0):
            rho2 = 1-(((P[2]*(1-P[3])).sum())/((P[1]*(1-P[3])).sum()))
        else:
            rho2 = 0.0

        if ((df[1]*P[1]*(1-P[2])*(1-P[3])).sum() > 0):
            rho3 = ((df[1]*P[1]*(1-P[2])*(1-P[3])).sum())/((P[1]*(1-P[2])*(1-P[3])).sum())/df[0].max()
        else:
            rho3 = 0.0

        if (df[2].sum()>0):
            rho4 = (df[2].sum())/(P[2].sum())/df[2].max()
        else:
            rho4=0.0

        rho5 = (stats.spearmanr(df[1],df[2])[0]/2)+0.5

        U2 = (df[2]==1).astype(int)
        U3 = (df[3]==0).astype(int)

        if ((U2*U3).sum()>0):
            rho6 = ((U2*U3).sum())/P[2].sum()
        else:
            rho6=0.0

        # rho7
        n_strings = 0
        s = 0
        for k in range(0, len(U2)):
            if (U2[k]*U3[k] == 0 and df_graph.has_node(k)):
                for l in list(nx.neighbors(df_graph,k)):
                    s += U2[l]*U3[l]
            #print(U2[i],U3[i])
        #lsprint(s)
        n_strings = s/2
        
        
        if (n_strings>0):
            rho7 = 1.0 - (n_strings/((U2*U3).sum()))
        else:
            rho7 = 1.0
        ###
        
        if (df[3].sum()>0):
            rho8 = df[3].sum()/(3*((df[2].sum())+(df[3].sum())))
        else:
            rho8 = 0.0

        if (P[3].sum()>0):
            rho9 = 1 - (P[3].sum()/df[3].sum())
        else:
            rho9=0.0

        if (P[3].sum()):
            rho10 = (P[3].sum())/(P[0].sum())
        else:
            rho10=0.0

        if ((df[0]*P[3]).sum()>0):
            rho11 = ((P[3]*U2).sum())/(P[3].sum())
        else:
            rho11=0.0

        if ((df[0]*P[3]).sum()>0):
            rho12 = ((df[0]*P[3]).sum())/(P[3].sum())/df[0].max()
        else:
            rho12=0.0
        
        rho = pd.DataFrame([[1,2,3,4,5,6,7,8,9,10,11,12],[rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8, rho9, rho10, rho11, rho12]]).T

        rho[0]=rho[0].astype(int)

        fullpath = path_to_outdir+f_network+'.orca.rho'

        rho.to_csv(fullpath, sep=' ',header=False, index=False)

def main(argv):
    path_to_orca = argv[1]
    path_to_outdir = argv[2]
    path_to_network = path_to_orca
    if len(argv) > 3:
        path_to_network = argv[3]

    #to the directory containing .ocount, -network.txt, .rho
    compute_rho(path_to_network,path_to_orca,path_to_outdir) 

if __name__ == "__main__":
    main(sys.argv)
