# Pathway Reconciliation

**Target Conference:** [Pacific Symposium on Biocomputing](https://psb.stanford.edu/)

**Target Session:** [Human Intrigue: Meta Analysis Approaches for Big Questions with Big Data](https://psb.humanintrigue.com/)

**Target Deadline:** August 2, 2021 ([Instructions for Authors](https://psb.stanford.edu/psb-online/psb-submit/))

**[Overleaf Document](https://www.overleaf.com/project/60edc899603043083d96cab6)**

**[Paper GitHub Repo](https://github.com/annaritz/2022-PSB-Human-Intrigue-Graphlets)**

Code repo ported from **[graphlet-tools](https://github.com/tobiasrubel/graphlet-tools)**

Some code related to parsing datasets ported from **[pathway-connectivity](https://github.com/annaritz/pathway-connectivity)** and **[pathway-parsers](https://github.com/Reed-CompBio/pathway-parsers)**.

## Directory Contents

- **src**: source code

- **networks**: parsed networks, null model networks, and interactomes.

- **graphlets**: orbit, graphlet, and rho values for networks.

- **regression**: regression files (ground truth and predictions).

- **out**: output visualization and analysis files.


## Reproducability

To reproduce plots in paper, as well as all attendant data, just run bash pipeline.sh. Run these commands first:

```
git submodule init
git submodule update
cd src/graphlets/ORCA/orca/
g++ -O2 -std=c++11 -o orca.exe orca.cpp
```

See the `README.md` file in `orca/` for more info.
