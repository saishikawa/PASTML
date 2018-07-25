# HIV-1A in Albania

This data set comes from the study by [Chevenet *et al.*, 2013](https://doi.org/10.1093/bioinformatics/btt010), which adopted the data from [Salemi *et al.*, 2008](https://doi.org/10.1371/journal.pone.0001390).

## Data

It contains a rooted tree with 152 tips (Albanian.tree.152tax.nwk), 
and a comma-separated annotation file describing the Country of sampling for each tree tip (state_Country.csv).

## Ancestral character reconstruction (ACR) analysis

To run the ACR of the Country with the state evolution model: F81, 
and the maximum likelihood prediction method: marginal approximation, run:
```bash
pastml -t Albania/data/pastml/Albanian.tree.152tax.nwk -a Albania/data/pastml/Country/F81/joint/state_Country.csv -m F81 -p marginal_approx
```

Other __models__ are also available: JC and EFT (which stands for "estimate frequencies from tips", not recommended),
as well as other __maximum likelihood prediction methods__ (marginal, max_posteriori, and joint), 
and __parsimony methods__ (downpass, acctran, or deltran). 
Parsimonious methods do not depend on branch lenghts so you do not need to specify an evolutionary model for them.