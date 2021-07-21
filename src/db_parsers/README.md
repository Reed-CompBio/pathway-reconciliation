# Pathway Database Parsers & Associated Code

This directory contains code to generate network files for signaling pathway databases. It also contains code to find corresponding pathways across databases.

## PathwayCommons

The first pathway databases were parsed from  [PathwayCommons](https://www.pathwaycommons.org/). These and all other downloaded files are in `infiles/`.

### Panther

Download `PathwayCommons11.panther.hgnc.txt.gz` from [PC2v11](https://www.pathwaycommons.org/archives/PC2/v11/) and unzip it.  

```
python3 parse_pathwaycommons.py -i infiles/PathwayCommons11.panther.hgnc.txt -o ../../networks/dbs/
```

### NCI-PID

Download `PathwayCommons11.pid.hgnc.txt.gz` from [PC2v11](https://www.pathwaycommons.org/archives/PC2/v11/) and unzip it.  

```
python3 parse_pathwaycommons.py -i infiles/PathwayCommons11.pid.hgnc.txt -o ../../networks/dbs/
```

### INOH

Download `PathwayCommons11.inoh.hgnc.txt.gz` from [PC2v11](https://www.pathwaycommons.org/archives/PC2/v11/) and unzip it.  

```
python3 parse_pathwaycommons.py -i infiles/PathwayCommons11.inoh.hgnc.txt -o ../../networks/dbs/

### PathBank?

## NDEx

The next pathway databases are from [NDEx](https://home.ndexbio.org/index/)

### SIGNOR

### CausalBioNet

## Original Pathway Databases

The final set of pathway databases are parsed from their original form. While they are available in some form from PathwayCommons, there are issues with each one (not all NetPath pathways area available, only KEGG metabolic pathways are availabile, and Reactome pathways are not parsed according to the hierarchy).

### NetPath

### KEGG

### Reactome
