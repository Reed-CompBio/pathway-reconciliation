# Pathway Database Parsers & Associated Code

This directory contains code to generate network files for signaling pathway databases. It also contains code to find corresponding pathways across databases.

All pathways must pass these minimal filters:
- All nodes must be proteins (e.g., have a UniProtKB identifier) unless they are "collapsed".
- There must be at least 10 interactions

## PathwayCommons

The first pathway databases were parsed from  [PathwayCommons](https://www.pathwaycommons.org/). These and all other downloaded files are in `infiles/`.

### Panther

Download `PathwayCommons11.panther.hgnc.txt.gz` from [PC2v11](https://www.pathwaycommons.org/archives/PC2/v11/) and unzip it in `infiles/`.

```
python3 parse_pathwaycommons.py -i infiles/PathwayCommons11.panther.hgnc.txt -o ../../networks/dbs/
```

30 too small, 0 nonhuman, 0 identical. 92 TOTAL parsed.

### NCI-PID

Download `PathwayCommons11.pid.hgnc.txt.gz` from [PC2v11](https://www.pathwaycommons.org/archives/PC2/v11/) and unzip it in `infiles/`.  

```
python3 parse_pathwaycommons.py -i infiles/PathwayCommons11.pid.hgnc.txt -o ../../networks/dbs/
```

10 too small, 0 nonhuman, 0 identical. 202 TOTAL parsed.

### INOH

Download `PathwayCommons11.inoh.hgnc.txt.gz` from [PC2v11](https://www.pathwaycommons.org/archives/PC2/v11/) and unzip it in `infiles/`. In addition, pathways labeled as Mammal, C. elegans, Drosophila, etc., are duplicate pathways and should be ignored.  37 pathways are ignored.

```
python3 parse_pathwaycommons.py -i infiles/PathwayCommons11.inoh.hgnc.txt -o ../../networks/dbs/
```

57 too small, 26 nonhuman, 37 identical. 113 TOTAL parsed.

## NDEx

The next pathway databases are from [NDEx](https://home.ndexbio.org/index/)

## Original Pathway Databases

The final set of pathway databases are parsed from their original form. While they are available in some form from PathwayCommons, there are issues with each one (not all NetPath pathways area available, only KEGG metabolic pathways are available, and Reactome pathways are not parsed according to the hierarchy).

### NetPath

Original 32 pathways from [NetPath](http://www.netpath.org/), used in previous publications such as [PathLinker](https://www.nature.com/articles/npjsba20162), are stored in `infiles/netpath`.  The script below makes the pathways undirected and two-column files to drop in `../../networks/dbs/`.

```
python3 parse_netpath.py
```

1 too small. 31 TOTAL parsed.

### KEGG

KEGG files are parsed according to an older [KEGG Pathway Parser](https://github.com/Reed-CompBio/pathway-parsers) which downloads KGML files and generates two types of files: `collapsed` files which includes protein complexes and protein families as individual nodes and `expanded` files that expand these entities into the individual proteins. These KGML files were downloaded Jul 14 2020 and moved to `infiles/kegg`.

The script below renames the pathways according to the `_pathway_list.txt` file in the KEGG directory, makes the edges undirected and two-columns. Drops them in `../../networks/kegg/collapsed/` and `../../networks/kegg/expanded`.

```
python3 parse_kegg.py
```

This file uses `mapping.txt` (downloaded Name & Uniprot from [HGNC mapper](https://www.genenames.org/download/custom/)). Also note that since we require at least 10 edges, this means that collapsed and expanded are slightly different in terms of the number of pathways. THis shouldn't affect the corresponding pathway analysis though.

- **collapsed:** 95 too small. 242 TOTAL parsed.
- **expanded:** 70 too small. 267 TOTAL parsed.

### Reactome

**TODO**

## Databases not Parsed
- PathBank
- SIGNOR
- CausalBioNet
- Reactome

# Getting Corresponding Pathways

```
python3 get_corresponding_pathways.py
```

This script generates, for every pathway, the best choice for some other pathway in every other database. It does so using the following rules:

## Candidate Pathways

- Two pathways are candidates if they share at least one term in the name, after removing domain specific general terms and stop words.
- Two pathways are candidates if they share at least one protein in common. _Note:_ `kegg_expanded` is used to compare proteins. Corresponding pathways are then automatically transferred to `kegg_collapsed` since they are derived from the same dataset.

For pathway _a_ from database _A_, the top corresponding pathway _b_ from database _B_ is a candidate pathway with the largest Jaccard overlap of the nodes (divided by the size of pathway _a_).

The file `jaccard.txt` records the top pick for every pathway _a_ from database _A_ compared to database _B_.  The file also contains the tokens in common for the names.  Note that this is assymetric.

## Corresponding Pathways

One can consider a graph, where the nodes are all pathways and there is a directed edge from pathway _a_ to pathway _b_ if _b_ is the top corresponding pathway of _a_ out of all pathways in database _B_.  However, this is a weakly connected graph with 553 edges.  

Instead, we will do something a little more stringent.  We define an undirected graph where there is an edge between pathway _a_ and pathway _b_ if _b_ is the top corresponding pathway of _a_ **AND** _a_ is the top corresponding pathway of _b_.  There are 142 such edges in this graph with 64 connected components. These are written in the file `conncomps.txt`.

Of the connected components, 23 connected components have members that are from at least 3 databases. We use these for our set of corresponding pathways.

## Manual Curation
