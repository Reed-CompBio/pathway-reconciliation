# AMI plots to compare real pathways with RW-Induced networks.
# argument is specified to tell the program whether we want to 
# compare corresponding pathways or for a specific database or for all
# Example run:
# python plot_ami_random.py corresponding
# or
# python plot_ami_random.py all 
# or database name like:
# python plot_ami_random.py netpath

import sys
sys.path.append("..")
import colors

import os
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_mutual_info_score
import glob


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
    #plt.yscale('log')
    plt.xlabel('# clusters')
    plt.ylabel('AMI')
    plt.savefig(svname)
    print(svname.split('/')[-1])

def plot_combined_ami(values,nets,outname):
    plt.clf()
    max_val = 0
    min_val = 1
    for k in values:
        max_val = max(max_val,max(values[k][1]))+0.01
        min_val = min(min_val,min(values[k][1]))-0.01
        max_x = max(values[k][0])
        min_x = min(values[k][0])

    plt.plot(values['gcount'][0],values['gcount'][1],'o-',color=colors.COLORS['darkblue'],ms=5,label='Graphlet Count')
    plt.plot(values['rho'][0],values['rho'][1],'s-',color=colors.COLORS['lightblue'],ms=5,label='GHuST Coefficient')
    
    plt.legend(loc=1,fontsize=8)
    plt.ylim([min_val,max_val])
    plt.grid()
    plt.xlabel('# clusters')
    plt.ylabel('AMI')
    if (nets!='all'):
        plt.title(nets)
    plt.tight_layout()
    plt.savefig(outname)
    print('Saved to %s' % (outname))

def aggregate_coeff(correspondence_file,path,ext,nets,net_type):

    cp = pd.read_csv(correspondence_file,sep='\t',index_col=0)
    dbs = list(cp.columns.values)
    pathways = list(cp.index.values)

    if (nets == 'corresponding'):

        if (net_type == 'original'):
            net = '-network.txt.orca.'
        elif (net_type == 'random-empirical'):
            #pick one of the random networks (0) for each original network
            net = '-network-RandomWalkerInduced-0-edges-network.txt.orca.'

        #elif (net_type == 'random-rewiring'):

        df=None
        for i in pathways:
            for j in dbs:
                if (pd.notnull(cp.loc[i,j])):
                    f = path+cp.loc[i,j][:-4]+net+ext
                    if (os.path.isfile(f)):
                        print(f)
                        df1=pd.read_csv(f,sep =' ', header=None)
                        df1=pd.DataFrame(df1[1])
                        df1=df1.rename(columns={1:j+'-'+i}).T
                        if df is not None:
                            df = df.append(df1)
                        else:
                            df = df1
        if df is not None:
            df['type']=net_type


    elif (nets == 'all'):
        dr = glob.glob(path+'/*/')
        dbs = [x.split('/')[-2] for x in dr]
        #print(dbs)
        df=None
        for j in dbs:

            if (net_type == 'original'):
                i = path+j+'/'
                filenames = glob.glob(i+'*-network.txt.orca.'+ext)
            elif (net_type == 'random-empirical'):
                i = path+j+'/'
                filenames = glob.glob(i+'*-network-RandomWalkerInduced-0-edges-network.txt.orca.'+ext)
                #pick one of the random networks (0) for each original network

            #elif (net_type == 'random-rewiring'):



            for f in filenames:

                if (os.path.isfile(f)):
                    print(f)
                    df1=pd.read_csv(f,sep =' ', header=None)
                    df1=pd.DataFrame(df1[1])
                    df1=df1.rename(columns={1:j+'-'+i}).T
                    if df is not None:
                        df = df.append(df1)
                    else:
                        df = df1
            if df is not None:
                df['type']=net_type

    else:
        dbs = [nets]
        df=None
        for j in dbs:

            if (net_type == 'original'):
                i = path+j+'/'
                filenames = glob.glob(i+'*-network.txt.orca.'+ext)
            elif (net_type == 'random-empirical'):
                i = path+j+'/'
                filenames = glob.glob(i+'*-network-RandomWalkerInduced-0-edges-network.txt.orca.'+ext)
                #pick one of the random networks (0) for each original network

            #elif (net_type == 'random-rewiring'):



            for f in filenames:

                if (os.path.isfile(f)):
                    print(f)
                    df1=pd.read_csv(f,sep =' ', header=None)
                    df1=pd.DataFrame(df1[1])
                    df1=df1.rename(columns={1:j+'-'+i}).T
                    if df is not None:
                        df = df.append(df1)
                    else:
                        df = df1
            if df is not None:
                df['type']=net_type

    return df


def main(argv):

    nets = argv[1] # corresponding or the particular database (e.g. netpath) or all
    correspondence_file = 'corresponding-top-picks-7pathways-6overlap-withcollapsed.txt'
    path_orig = r'../../graphlets/dbs/'
    path_cp = r'../../networks/dbs/'
    path_rand = r'../../graphlets/null-models/random-empirical/'
    path_out = r'../../out/'
    correspondence_file=path_cp+correspondence_file

    values = {}
    for ext in list(['rho','gcount']):
        values[ext] = {}
        #original
        net_type = 'original'

        df = aggregate_coeff(correspondence_file,path_orig,ext,nets,net_type)

        #random
        net_type = 'random-empirical'

        df1 = aggregate_coeff(correspondence_file,path_rand,ext,nets,net_type)

        df = df.append(df1)
        ground_truth=df['type']
        df = df[df.columns[0:len(df.columns)-1]]

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
        amiplnm = path_out+'AMI_RW-Induced-'+ext+'-'+nets+'.pdf' 
        values[ext] = (irange,amirange)
        plotami(irange,amirange,amiplnm)

    plot_combined_ami(values, nets, '%s/Combined_AMI_RW-Induced-%s.pdf' % (path_out, nets))



if __name__ == "__main__":
    main(sys.argv)
