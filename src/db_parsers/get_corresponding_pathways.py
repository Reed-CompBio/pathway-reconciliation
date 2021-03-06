## gets corresponding pathways.
import glob
import re
import sys
import networkx as nx


domain_specific = set(['signaling','signalling','pathway','receptor','events','heterodimer','dimer','activation','inhibition','negative','positive','diagram','regulation','network','transcription','effects','factor','cascade','cancer','disease','action','acid','hormone','growth','guidance','cell','cells','cellular','biosynthesis','complement','contraction','complex','cycle','degredation','dna','downstream','effectors','expression','exchange','family','gene','interaction','interactions','kinase','phosphatase','longterm','migration','multiple','nuclear','outer','phosphate','phosphorylation','posttranslational','presentation','processing','protein','regulate','regulates','regulated','response','role','roles','signal','signals','species','specific','stabilization','stability','surface','target','targets','channels','transduction','transport','proteins'])
stopwords = set(['of','by','and','the','on','from','to','in','other','with','through','non','a','via'])

PATHWAY_FILES = {'netpath':'../../networks/dbs/netpath/pathway-names.txt',
    'kegg_expanded':'../../networks/dbs/kegg_expanded/pathway-names.txt',
    'panther':'../../networks/dbs/panther/pathway-names.txt',
    'inoh':'../../networks/dbs/inoh/pathway-names.txt',
    'pid':'../../networks/dbs/pid/pathway-names.txt',
    'pathbank':'../../networks/dbs/pathbank/pathway-names.txt',
    'signor_expanded':'../../networks/dbs/signor_expanded/pathway-names.txt'}
    #'reactome':'../../networks/dbs/reactome/pathway-names.txt'}

DB_ORDER = ['netpath','pid','panther','inoh', 'pathbank', 'signor_expanded','kegg_expanded']

def main():
    pw_order,names,nodes = read_info()
    edges,top_picks = get_top_picks(pw_order,names,nodes)
    corresponding_pathways(edges,top_picks)

def read_info():
    pw_order = {}
    nodes = {}
    names = {}
    for db in DB_ORDER:
        prefix = PATHWAY_FILES[db].replace('pathway-names.txt','')
        pw_order[db] = []
        nodes[db] = {}
        names[db] = {}
        with open(PATHWAY_FILES[db]) as fin:
            for line in fin:
                row = line.strip().split('\t')
                #print(row)
                pw_name = row[0]
                if skip(row[2]):
                    continue
                pw_order[db].append(pw_name)
                names[db][pw_name] = tokenize(row[2])
                nodes[db][pw_name] = set()
                with open(prefix+row[0]) as fin2:
                    for line2 in fin2:
                        row2 = line2.strip().split()
                        nodes[db][pw_name].add(row2[0])
                        nodes[db][pw_name].add(row2[1])
    return pw_order,names,nodes

def get_top_picks(pw_order,names,nodes):
    edges = set()
    out = open('jaccard-%dpathways-%doverlap.txt' % (len(DB_ORDER),THRES),'w')
    out.write('#db1\tpw1\tdb2\tpw2\tname_jaccard\tnode_jaccard\n')
    top_picks = {}
    for i in range(len(DB_ORDER)):
        dbi = DB_ORDER[i]
        top_picks[dbi] = {}
        for j in range(len(DB_ORDER)):
            if i == j:
                continue
            dbj = DB_ORDER[j]
            top_picks[dbi][dbj] = {}
            print(dbi,dbj)
            for pi in pw_order[dbi]:
                best_pw = None
                best_pw_jaccard = None
                best_pw_name = None
                for pj in pw_order[dbj]:
                    jaccard_name = names[dbi][pi].intersection(names[dbj][pj])
                    jaccard_nodes = jaccard(nodes[dbi][pi],nodes[dbj][pj])
                    if len(jaccard_name) > 0 and jaccard_nodes > 0:
                        if best_pw == None or jaccard_nodes > best_pw_jaccard:
                            best_pw_jaccard = jaccard_nodes
                            best_pw_name = jaccard_name
                            best_pw = pj
                if best_pw != None:
                    print('-->',pi,'  |  ',best_pw,'  |  Name Jaccard:',best_pw_name,'Node Jaccard:',best_pw_jaccard)
                    out.write('%s\t%s\t%s\t%s\t%s\t%f\n' % (dbi,pi,dbj,best_pw,'|'.join([w for w in best_pw_name]),best_pw_jaccard))
                    #print('    adding (%s_%s,%s_%s)' % (dbi,pi,dbj,best_pw))
                    edges.add(('%s_%s' % (dbi,pi),'%s_%s' % (dbj,best_pw)))
                    top_picks[dbi][dbj][pi]={'pij':best_pw,'name':best_pw_name,'jaccard':best_pw_jaccard}

    out.close()
    print('wrote to jaccard')
    return edges,top_picks

def corresponding_pathways(edges,top_picks):
    # Step 1. Build graph from symmetric top picks
    bidirected_edges = set()
    for u,v in edges:
        if (v,u) in edges and (v,u) not in bidirected_edges:
            bidirected_edges.add((u,v))

    print(len(edges),len(bidirected_edges))
    G = nx.from_edgelist(bidirected_edges,create_using=nx.Graph)
    print('%d nodes and %d edges.' % (G.number_of_nodes(),G.number_of_edges()))

    ## Step 2. Get connected components
    conncomps = []
    out = open('conncomps-%dpathways-%doverlap.txt'  % (len(DB_ORDER),THRES),'w')
    for conncomp in nx.connected_components(G):
        out.write('%d\t%s\n' % (len(conncomp),'\t'.join([c for c in conncomp])))
        dbs = set([c.split('_')[0] for c in conncomp])
        if len(dbs) >= THRES :
            conncomps.append([c for c in conncomp])
    out.close()
    print('wrote to conncomps')
    print('%d conncomps have %d or more symmetric top pick members' % (len(conncomps),THRES))

    ## make corresponding-top-picks.txt file
    outfile = '../../networks/dbs/corresponding-top-picks-%dpathways-%doverlap-ORIG.txt' % (len(DB_ORDER),THRES)
    out = open(outfile,'w')
    out.write('#pathway\t'+'\t'.join(DB_ORDER)+'\n')
    for conncomp in sorted(conncomps,key=len):
        row = {DB:[] for DB in DB_ORDER}
        row['kegg_collapsed'] = []
        for c in conncomp:
            if 'kegg' in c or 'signor' in c:
                this_db = '_'.join(c.split('_')[:2])
                this_pw = '_'.join(c.split('_')[2:])
            else:
                this_db = c.split('_')[0]
                this_pw = '_'.join(c.split('_')[1:])
            row[this_db].append('%s/%s' % (this_db,this_pw))
        out.write('NAME')
        for db in DB_ORDER:
            if row[db] == []:
                out.write('\tNaN')
            else:
                out.write('\t'+'|'.join(row[db]))
        out.write('\n')
    out.close()
    print('wrote to',outfile)
    return

def skip(text):
    text = text.lower()
    text = re.sub(r'[^\w\s]', '',text) # remove punctuation
    to_ignore = set(['cancer','carcinoma','disease','inflammation','viral','virus','carcinogenesis','diabetes','asthma','measles','resistance','infection','xenopus','mouse','elegans','drosophila','mammal','metabolism','metabolic','syndrome'])
    for word in to_ignore:
        if word in text:
            return True
    return False

def tokenize(text):
    # tips from https://dylancastillo.co/nlp-snippets-clean-and-tokenize-text-with-python/
    text = text.lower()
    text = re.sub(r'[^\w\s]', '',text) # remove punctuation
    text = text.replace(' cell ','_cell ') # for B_cell, T_cell, etc.
    text = set([t for t in text.split() if len(t) > 1]) # ignore single-letter tokens
    text = text.difference(domain_specific)
    text = text.difference(stopwords)
    #print(text)
    return text

def jaccard(set1,set2):
    if len(set1)==0:
        return 0
    return len(set1.intersection(set2))/len(set1)

if __name__ == '__main__':
    global THRES
    if len(sys.argv)==1:
        THRES = 4
    else:
        THRES = int(sys.argv[1])
    print('threshold is',THRES)
    main()
