import glob
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot


PATHWAY_FILES = {
    'netpath':'../../networks/dbs/netpath/*',
    'panther':'../../networks/dbs/panther/*',
    'inoh':'../../networks/dbs/inoh/*',
    'pid':'../../networks/dbs/pid/*',
    'kegg_collapsed':'../../networks/dbs/kegg_collapsed/*',
    'kegg_expanded':'../../networks/dbs/kegg_expanded/*',
    'reactome':'../../networks/dbs/reactome/*'}

lens = {}
fig = plt.figure(figsize=(4,8))
max_len = 0

i = 1
for p in PATHWAY_FILES:

    files = glob.glob(PATHWAY_FILES[p])
    lens[p] = []
    for f in files:
        if 'pathway-name' in f:
            continue
        # from https://www.oreilly.com/library/view/python-cookbook/0596001673/ch04s07.html
        lens[p].append(len(open(f).readlines()))
    max_len = max(max_len,max(lens[p]))

    print(p,len(lens[p]))
    plt.subplot(len(PATHWAY_FILES),1,i)
    plt.title('%s (n=%d, max_len=%d)' % (p,len(lens[p]),max(lens[p])))
    pyplot.hist(lens[p],np.linspace(0,max(1000,max(lens[p])),100),alpha=.75,edgecolor='k')
    if p == 'reactome':
        pyplot.plot([1000,1000],[0,500],'-r')
    else:
        pyplot.plot([1000,1000],[0,10],'-r')

    i+=1

plt.tight_layout()
plt.savefig('pathway_sizes.png')
print('wrote to pathway_sizes.png')
