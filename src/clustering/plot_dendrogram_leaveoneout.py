from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import fcluster
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.patches as patches
import random
sys.path.append('..')
import colors

def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)



def get_linkage(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)
    return linkage_matrix
    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)


def cluster__(df,dest,exclude_txt,color_mapping):
    global L
    X = df.T.values
    print(X)
    link = 'average'
    affin = 'cosine'
    #fig = plt.figure(1, figsize=(14, 2))
    fig = plt.figure(1, figsize=(14, 2))
    ax = fig.add_subplot(111)
    model = AgglomerativeClustering(distance_threshold=0, n_clusters=None,affinity=affin,linkage=link)
    model = model.fit(X)
    #lbls = ['-'.join(x.split('-')[:-1]) for x in df.columns]
    lbls = [x for x in df.columns]
    print(lbls)
    plot_dendrogram(model, labels=lbls,leaf_rotation=90, leaf_font_size=8,color_threshold=0,above_threshold_color='k')#truncate_mode='lastp', p=k,)#leaf_label_func=llf)
    #add colors
    _,l = plt.xticks()
    l = [x.get_text() for x in l]
    x = 0
    #cs = ['#20A4F3','#E84A7F','#FF8A5B','#2E294E','#7FB685','#FFFD82','#E0A890','#C7DFC5','#C1DBE3']
    #cs = list(colors.COLORS.values())
    #cs = ['#20A4F3','#AAAAAA']
    print('color_mapping:',color_mapping)
    num_to_color = -1
    for c in color_mapping:
        if exclude_txt.split('_')[0]==c:
            num_to_color = color_mapping[c]
    num_to_color=str(num_to_color)
    print('num to color:',num_to_color)
    num = 0 # num colored
    tot = 0
    x_ticks = []
    x_labels = []
    for i in l:
        if i.split('_')[-1] == num_to_color:
            c = colors.COLORS['darkblue']
            print(i)
            num+=1
            x_ticks.append(x+5)
            x_labels.append('\n'.join(i.split('_')[1:-1]))
        else:
            c = colors.COLORS['gray']
        tot+=1
        ax.add_patch(
        patches.Rectangle(
            xy=(x,0),  # point of origin.
            width=10,
            height=1,
            linewidth=0,
            color=c,
            alpha=0.8,
            fill=True))
        x+=10
    print('%d total colored' % (num))
    #plt.ylim(0.6, 1)
    plt.xticks(x_ticks,x_labels)
    #plt.xticks(fontsize=16)
    #plt.yticks(fontsize=16)

    #plt.xticks([], [])
    #plt.xticks(ticks=[5+x*10 for x in range(len(df.columns))],labels=['-'.join(x.split('-')[:-1]) for x in df.columns])
    #plt.tight_layout()
    #plt.axis('off')
    plt.title(exclude_txt[:-1])
    plt.tight_layout()

    plt.savefig(dest,bbox_inches='tight')
    #plt.show()

def read_vecs(path: str) -> pd.DataFrame:
    """
    :path    to csv
    :returns DataFrame of column vectors
    """
    df = pd.read_csv(path,index_col=0)#.apply(lambda x: 1-x)
    return df

def write__(classes,dest):
    for c in classes:
        with open(os.path.join(dest,'{}.txt'.format(c)),'w') as f:
            f.write('\n'.join(classes[c]))

def combine(lat):
    """
    :lat     list of paths to gcount files
    :returns df of all graphlets
    """
    return pd.concat([pd.read_csv(x,index_col=0,sep='\s+',names=[x.split('/')[-1].split('.')[0]+'_'+x.split('/')[-2]]) for x in lat], axis=1, sort=True)

def combine2(lat,pw):
    """
    :lat     list of paths to gcount files
    :returns df of all graphlets
    """
    return pd.concat([pd.read_csv(x,index_col=0,sep='\s+',names=[pw+'_'+x.split('/')[-2]]) for x in lat], axis=1, sort=True)

def merge_and_label(lat):
    print('BEFORE:',[lat[i].columns for i in range(len(lat))])
    for i in range(len(lat)):
        lat[i].columns = ['{}_{}'.format(x,i) for x in lat[i].columns]
    print('AFTER:',[lat[i].columns for i in range(len(lat))])
    ndf = pd.concat(lat,axis=1)
    return ndf

def process_inputs(lat,pathwayassoc,color_by_pathway,exclude_txt) -> pd.DataFrame:
    """
    :lat              list of gcount files
    :pathwayassoc     association between networks and databases/pathways
    :color_by_pathway whether to color by pathway or by db
    :returns          labeled dataframe for plotting
    """
    print('\nprocessing inputs...')
    pathwaykey = pd.read_csv(pathwayassoc,index_col=0,sep='\t')
    #lets filter out the elements of the lat which are not in the pathwaykey. Uncomment this to see everything vs everything.
    print(lat[0])
    print('/'.join(lat[0].split('/')[-2:]))
    print(pathwaykey)

    #return
    lat = [x for x in lat if '-'.join('/'.join(x.split('/')[-2:]).split('-')[:-1]).replace(exclude_txt,'') in list(pathwaykey.applymap(lambda x: x.split('.')[0] if isinstance(x,str) else np.nan).stack())]
    print(lat)
    color_mapping = {}
    if color_by_pathway:
        print(pathwaykey)
        storage = []
        n = 0
        for pw in pathwaykey.index:
            print(pw)
            if pw == exclude_txt.split('_')[0]:
                num_to_exclude=n
            color_mapping[pw]=n
            print(pw)
            if pw == 'p38/MAPK':
                color_mapping['p38'] = n
            if pw == 'TNFalpha/Fas':
                color_mapping['TNFalpha'] = n
            n+=1
            corresponding = [pathwaykey.loc[pw][db] for db in pathwaykey.columns]
            corresponding = [a.replace('/','/'+exclude_txt).replace('.txt','-network.txt') for a in corresponding if isinstance(a,str)]
            pws = [x for x in lat if any([corr in x for corr in corresponding])]
            if pws != []:
                dbdf = combine2(pws,pw)
                storage.append(dbdf)
        final_df = merge_and_label(storage)
    else:
        storage = []
        for db in pathwaykey.columns:
            dbpws = [x for x in lat if x.split('/')[-2] == db]
            print('associating {} with {}'.format(db,dbpws))
            if dbpws != []:
                dbdf = combine(dbpws)
                storage.append(dbdf)
        final_df = merge_and_label(storage)
    return final_df, color_mapping




def main(argv):
    """
    :data        list of gcount files
    :exclude_txt  text to exclude
    :out         path and filename to save the dend as
    :side-effect plots a dendrogram of the data

    """
    pathwayassoc = '../../networks/dbs/corresponding-top-picks-7pathways-4overlap-withcollapsed.txt'
    color_by_pathway = True
    data = argv[1:-2]
    data = [x for x in data if not 'pathway-names.txt' in x]
    print(data)
    processed_data,color_mapping = process_inputs(data,pathwayassoc,color_by_pathway,argv[-2])
    print(processed_data)
    out = argv[-1]
    clusters = cluster__(processed_data,out,argv[-2],color_mapping)

if __name__ == "__main__":
    main(sys.argv)
