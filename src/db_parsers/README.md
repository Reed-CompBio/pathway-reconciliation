# Pathway Database Parsers & Associated Code

This directory contains code to generate network files for signaling pathway databases. It also contains code to find corresponding pathways across databases.

## PathwayCommons

The first pathway databases were parsed from  [PathwayCommons](https://www.pathwaycommons.org/). These and all other downloaded files are in `infiles/`.

### Panther

Download `PathwayCommons11.panther.hgnc.txt.gz` from [PC2v11](https://www.pathwaycommons.org/archives/PC2/v11/) and unzip it in `infiles/`.

```
python3 parse_pathwaycommons.py -i infiles/PathwayCommons11.panther.hgnc.txt -o ../../networks/dbs/
```

### NCI-PID

Download `PathwayCommons11.pid.hgnc.txt.gz` from [PC2v11](https://www.pathwaycommons.org/archives/PC2/v11/) and unzip it in `infiles/`.  

```
python3 parse_pathwaycommons.py -i infiles/PathwayCommons11.pid.hgnc.txt -o ../../networks/dbs/
```

### INOH

Download `PathwayCommons11.inoh.hgnc.txt.gz` from [PC2v11](https://www.pathwaycommons.org/archives/PC2/v11/) and unzip it in `infiles/`.

```
python3 parse_pathwaycommons.py -i infiles/PathwayCommons11.inoh.hgnc.txt -o ../../networks/dbs/
```

### PathBank

**TODO**

## NDEx

The next pathway databases are from [NDEx](https://home.ndexbio.org/index/)

### SIGNOR

**TODO**

### CausalBioNet

**TODO**

## Original Pathway Databases

The final set of pathway databases are parsed from their original form. While they are available in some form from PathwayCommons, there are issues with each one (not all NetPath pathways area available, only KEGG metabolic pathways are available, and Reactome pathways are not parsed according to the hierarchy).

### NetPath

Original 32 pathways from [NetPath](http://www.netpath.org/), used in previous publications such as [PathLinker](https://www.nature.com/articles/npjsba20162), are stored in `infiles/netpath`.  The script below makes the pathways undirected and two-column files to drop in `../../networks/dbs/`.

```
python3 parse_netpath.py
```

### KEGG

KEGG files are parsed according to an older [KEGG Pathway Parser](https://github.com/Reed-CompBio/pathway-parsers) which downloads KGML files and generates two types of files: `collapsed` files which includes protein complexes and protein families as individual nodes and `expanded` files that expand these entities into the individual proteins. These KGML files were downloaded Jul 14 2020 and moved to `infiles/kegg`.

The script below renames the pathways according to the `_pathway_list.txt` file in the KEGG directory, makes the edges undirected and two-columns. Drops them in `../../networks/kegg/collapsed/` and `../../networks/kegg/expanded`.

```
python3 parse_kegg.py
```

This file uses `mapping.txt` (downloaded Name & Uniprot from [HGNC mapper](https://www.genenames.org/download/custom/)).

### Reactome

**TODO**
