# Pathway Reconciliation

This repo contains code for the following paper:

**Reconciling Signaling Pathway Databases with Network Topologies.**
Tobias Rubel, Pramesh Singh, and Anna Ritz.
Pacific Symposium on Biocomputing (PSB) 2022. [PSB Proceedings Link.](https://psb.stanford.edu/psb-online/proceedings/psb22/rubel.pdf)

## Directory Contents

- **src**: source code

- **networks**: parsed networks, null model networks, and interactomes.

- **graphlets**: orbit, graphlet, and rho values for networks.

- **regression**: regression files (ground truth and predictions).

- **out**: output visualization and analysis files.


## Reproducibility

To reproduce plots in paper, as well as all attendant data, just run bash pipeline.sh. Run these commands first:

```
git submodule init
git submodule update
cd src/graphlets/ORCA/orca/
g++ -O2 -std=c++11 -o orca.exe orca.cpp
```

See the `README.md` file in `orca/` for more info.


### Developer Notes

Code repo ported from **[graphlet-tools](https://github.com/tobiasrubel/graphlet-tools)**

Some code related to parsing datasets ported from **[pathway-connectivity](https://github.com/annaritz/pathway-connectivity)** and **[pathway-parsers](https://github.com/Reed-CompBio/pathway-parsers)**.
