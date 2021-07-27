import glob
import numpy as np

DB_ORDER = ['netpath','pid','panther','inoh', 'pathbank', 'signor_expanded','kegg_expanded','signor_collapsed','kegg_collapsed']
DB_SOURCE = {'netpath':'NetPath~\\cite{kandasamy2010netpath}','pid':'PathwayCommons~\\cite{rodchenkov2020pathway}','panther':'PathwayCommons~\\cite{rodchenkov2020pathway}','inoh':'PathwayCommons~\\cite{rodchenkov2020pathway}', 'pathbank':'Pathbank~\\cite{wishart2020pathbank}', 'signor_expanded':'NDEx~\\cite{pillich2017ndex}','kegg_expanded':'KEGG~\\cite{kanehisa2021kegg}','signor_collapsed':'NDEx~\\cite{pillich2017ndex}','kegg_collapsed':'KEGG~\\cite{kanehisa2021kegg}'}
for key in DB_SOURCE:
    DB_SOURCE[key] = '\\small{%s}' % (DB_SOURCE[key])
DB_NAMES = {'netpath':'NetPath~\\cite{kandasamy2010netpath}','pid':'PID~\\cite{ncipid}','panther':'Panther~\\cite{mi2021panther}','inoh':'INOH\\cite{yamamoto2011inoh}', 'pathbank':'PathBank~\\cite{wishart2020pathbank}', 'signor_expanded':'\\textit{SIGNOR-expanded}~\\cite{licata2020signor}','kegg_expanded':'\\textit{KEGG-expanded}~\\cite{kanehisa2021kegg}','signor_collapsed':'SIGNOR~\\cite{licata2020signor}','kegg_collapsed':'KEGG~\\cite{kanehisa2021kegg}'}
DB_SCOPE = {'netpath':'Immune \& Cancer','pid':'Cancer','panther':'Primary Signaling','inoh':'Hierarchical Model', 'pathbank':'Model Organisms', 'signor_expanded':'Binary Causal','kegg_expanded':'Broad Focus','signor_collapsed':'Binary Causal','kegg_collapsed':'Broad Focus'}
DB_ORDER = sorted(DB_ORDER)


print('\\begin{table}[h]')
print('\\begin{tabular}{|lll|ccc|}\\hline')

#print('& & \\textbf{Source for} & &  & \\\\')
print(' \\textbf{Database} & \\textbf{Focus} & \\textbf{Parse Source} & \\textbf{$n$} & \\textbf{\\# Nodes} & \\textbf{\\# Edges} \\\\ \\hline')
for db in DB_ORDER:
    files = glob.glob('../../networks/dbs/%s/*.txt' % (db))
    n = 0
    edges = []
    nodes = []
    for pw in files:
        if 'pathway-name' in files:
            continue
        n+=1
        with open(pw) as fin:
            these_nodes = set()
            num_edges = 0
            for line in fin:
                num_edges+=1
                row = line.strip().split('\t')
                these_nodes.add(row[0])
                these_nodes.add(row[1])
            edges.append(num_edges)
            nodes.append(len(these_nodes))
    print('%s & %s & %s & $%d$ & $%d \pm %.2f$ & $%d \pm %.2f$ \\\\' % (DB_NAMES[db],DB_SCOPE[db],DB_SOURCE[db],n,np.mean(nodes),np.std(nodes),np.mean(edges),np.std(edges)))
    #print(DB_NAMES[db],DB_SOURCE[db],n,np.mean(nodes),np.std(nodes),np.mean(edges),np.std(edges))

these_nodes = set()
num_edges = 0
with open('../../networks/interactomes/All_Pathway_Commons.txt') as fin:
    for line in fin:
        num_edges+=1
        row = line.strip().split('\t')
        these_nodes.add(row[0])
        these_nodes.add(row[1])
    print('Interactome & Broad Focus & PathwayCommons~\\cite{rodchenkov2020pathway} & $1$ & $%d$ & $%d$ \\\\' % (len(these_nodes),num_edges))
