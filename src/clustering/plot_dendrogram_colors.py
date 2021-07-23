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
    fig = plt.figure(1, figsize=(14, 2))
    ax = fig.add_subplot(111)
    model = AgglomerativeClustering(distance_threshold=0, n_clusters=None,affinity=affin,linkage=link)
    model = model.fit(X)
    lbls = ['-'.join(x.split('-')[:-1]) for x in df.columns]
    plot_dendrogram(model, labels=lbls,leaf_rotation=90, leaf_font_size=8,color_threshold=0,above_threshold_color='k')#truncate_mode='lastp', p=k,)#leaf_label_func=llf)
    #add colors
    _,l = plt.xticks()
    l = [x.get_text() for x in l]
    x = 0
    c1='#20A4F3'
    c2='#E84A7F'
    for i in l:
        if 'Rand' in i:
            c = c1
        else:
            c = c2
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
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.xticks([], [])
    plt.tight_layout()
    #plt.axis('off')
    #plt.title('dendrogram of cos sim') 
    plt.savefig(dest)
   
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

def main(argv):
    """
    data = 'Test-Data/test-100.csv'
    clusters = cluster(read_vecs(data),90,True)
    write(clusters,'Test-Out')
    """
    data = argv[1]
    out = argv[2]
    clusters = cluster__(read_vecs(data),out,)

if __name__ == "__main__":
    main(sys.argv)

