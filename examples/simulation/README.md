# Simulation datasets
  
The dataset is used in the simulation study performed by Ishikawa et al., submitted to *Mol. Biol. Evol.*
  
## Data contained in each repository
  
__trees_1000_original__ : the pure-birth 1000-tips trees generated for 0.1 - 10 speciation substitution rate ratios  
__trees_1000_flat__ : same trees with those in 'trees_1000_original' except that their branch lengths are averaged (flattened)  
__Annotations_4__ : tip annotations (DNA-like data) generated on the original 1000-tips trees with the HKY model  
__Annotations_20__ : tip annotations (protein-like data) generated on the original 1000-tips trees with the JTT model  
__scenarios_4__ : the list of true ancestral states (A, T, G or C) for all internal nodes of the original tree, which were observed in the data generation  
__scenarios_20__ : the list of true ancestral states (20 amino-acids) for all internal nodes of the original tree, which were observed in the data generation  
__PASTML-v0.6.6.1__ : PASTML source code and binary to be used in the ancestral state reconstruction  
__Seq-Gen.v.1.3.3__ : Seq-Gen source code and binary to be used in the data simulation  
__Results_simulation__: output directories to contain the summary statistics for the accuracy of the ML- and Parsimony-based ancestral reconstruction. R commands are included to visualize the results.  
  
## To automatically analyze all data
  
```bash
./auto_analysis.sh
```
  
See our manuscript for details of the simulation procedure. A preprint will be available from bioRxiv soon.
