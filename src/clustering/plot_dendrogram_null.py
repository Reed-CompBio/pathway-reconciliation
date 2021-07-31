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
from colors import COLORS



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


def cluster__(df,dest):
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
    plot_dendrogram(model, labels=lbls,leaf_rotation=90, leaf_font_size=8,color_threshold=0,above_threshold_color='k')#truncate_mode='lastp', p=k,)#leaf_label_func=llf)
    #add colors
    _,l = plt.xticks()
    l = [x.get_text() for x in l]
    x = 0
    #cs = ['#20A4F3','#E84A7F','#FF8A5B','#2E294E','#7FB685','#FFFD82','#E0A890','#C7DFC5','#C1DBE3']
    cs = [COLORS['oceangray'],COLORS['orange']]
    for i in l:
        print(i)
        print(i.split('_')[-1])
        print('{} belongs with {}'.format(i,cs[int(i.split('_')[-1])]))
        c = cs[int(i.split('_')[-1])]
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

    #plt.ylim(0.6, 1)
    #plt.xticks(fontsize=16)
    #plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.xticks([], [])
    #plt.xticks(ticks=[5+x*10 for x in range(len(df.columns))],labels=['-'.join(x.split('-')[:-1]) for x in df.columns])
    #plt.tight_layout()
    #plt.axis('off')
    #plt.title('dendrogram of cos sim') 
    plt.savefig(dest)
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
    return pd.concat([pd.read_csv(x,index_col=0,sep='\s+',names=[x.split('/')[-1].split('.')[0]]) for x in lat], axis=1, sort=True)

def merge_and_label(lat):
    for i in range(len(lat)):
        lat[i].columns = ['{}_{}'.format(x,i) for x in lat[i].columns]
    ndf = pd.concat(lat,axis=1)
    return ndf

def process_inputs(lat,) -> pd.DataFrame:
    """
    :lat              list of gcount files
    :pathwayassoc     association between networks and databases/pathways
    :color_by_pathway whether to color by pathway or by db 
    :returns          labeled dataframe for plotting
    """
    #lets filter out the elements of the lat which are not in the pathwaykey. Uncomment this to see everything vs everything.
    random_nets = [x for x in lat if 'Rand' in x]
    pways = [x for x in lat if not x in random_nets]
    rdf = combine(random_nets)
    pdf = combine(pways)
    return merge_and_label([rdf,pdf])
    


def main(argv):
    """
    :data        list of gcount files
    :out         path and filename to save the dend as 
    :side-effect plots a dendrogram of the data 

    """
    data = argv[1:-1]
    data = [x for x in data if not 'pathway-names.txt' in x]
    print(data)
    processed_data = process_inputs(data)
    print(processed_data)
    out = argv[-1]
    clusters = cluster__(processed_data,out,)

if __name__ == "__main__":
    main(sys.argv)

