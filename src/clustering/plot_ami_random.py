# 
# Example run:
# python plot_ami.py corresponding-top-picks-7pathways-6overlap-withcollapsed.txt path/to/output/
# 

import os
import pandas as pd
import sys
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_mutual_info_score
#import glob


def cluster(df,numclusters):
    global model
    X = df.T.values
    link = 'average'
    affin = 'cosine'
    model = AgglomerativeClustering(distance_threshold=None, n_clusters=numclusters,affinity=affin,linkage=link)
    model = model.fit(X)
    mapping = dict(zip(df.columns,model.labels_))
    return model.labels_
    #lets reverse the mapping into a classification
    reversed_dict = dict()
    for k,v in mapping.items():
        try:
            reversed_dict[v].append(k)
        except KeyError:
            reversed_dict[v]=[k]
    return reversed_dict

def plotami(x,y,svname):
    plt.clf()
    plt.scatter(x,y)
    plt.grid()
    plt.xlabel('# clusters')
    plt.ylabel('AMI')
    plt.savefig(svname)
    print(svname.split('/')[-1])

def aggregate_coeff(correspondence_file,path,ext,net_type):

    if (net_type == 'original'):
        net = '-network.txt.orca.'
    elif (net_type == 'random-empirical'):
        #pick one of the random networks (0) for each original network
        net = '.txt.id-RandomWalkerInduced-0-edges-network.txt.orca.'

    #elif (net_type == 'random-rewiring'):


    cp = pd.read_csv(correspondence_file,sep='\t',index_col=0)
    dbs = list(cp.columns.values)
    pathways = list(cp.index.values)

    df=None
    for i in pathways:
        for j in dbs:
            if (pd.notnull(cp.loc[i,j])):
                f = path+cp.loc[i,j][:-4]+net+ext
                if (os.path.isfile(f)):

                    df1=pd.read_csv(f,sep =' ', header=None)
                    df1=pd.DataFrame(df1[1])
                    df1=df1.rename(columns={1:j+'-'+i}).T
                    df1['db']=j
                    df1['pathway']=i
                    if df is not None:
                        df = df.append(df1)
                    else:
                        df = df1
        if df is not None:
            df['type']=net_type
    return df


def main(argv):


    correspondence_file = argv[1]
    path_orig = r'../../graphlets/dbs/'
    path_cp = r'../../networks/dbs/'
    path_rand = r'../../graphlets/null-models/random-empirical/'
    path_out = r'../../out/'
    correspondence_file=path_cp+correspondence_file


    for ext in list(['rho','gcount']):
        #nets
        net_type = 'original'

        df = aggregate_coeff(correspondence_file,path_orig,ext,net_type)

        #random nets
        net_type = 'random-empirical'

        df1 = aggregate_coeff(correspondence_file,path_rand,ext,net_type)

        df = df.append(df1)
        ground_truth=df['type']
        df = df[df.columns[0:len(df.columns)-3]]

        if (ext=='gcount'):
            for i in range(0,len(df.index)):
                s=df.iloc[i,0]
                if s > 0:
                    df.iloc[i,0]=df.iloc[i,0]/s

                s=df.iloc[i,1:3].sum()
                if s > 0:
                    df.iloc[i,1:3]=df.iloc[i,1:3]/s

                s=df.iloc[i,3:9].sum()
                if s > 0:
                    df.iloc[i,3:9]=df.iloc[i,3:9]/s

                s=df.iloc[i,9:30].sum()
                if s > 0:
                    df.iloc[i,9:30]=df.iloc[i,9:30]/s


        df=df.T

        mn_nm = 2
        mx_nm = len(df.columns)
        irange = list(range(mn_nm,mx_nm+1))

        amirange = [adjusted_mutual_info_score(cluster(df,i),ground_truth) for i in irange]
        amiplnm = path_out+'by-network-type'+'-ami-'+ext+'.pdf' 
        plotami(irange,amirange,amiplnm)



if __name__ == "__main__":
    main(sys.argv)
