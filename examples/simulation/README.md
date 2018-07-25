# Simulation datasets and scripts
  
The dataset were used in the simulation study performed by Ishikawa et al., submitted to *Mol. Biol. Evol.*
  
## Data
  
__trees_1000_original__ : the pure-birth 1000-tips trees generated for 0.1 - 10 speciation substitution rate ratios  
__trees_1000_flat__ : same trees with *trees_1000_original* except that their branch lengths are averaged (flattened)  
__Annotations_4__ : tip annotations (DNA-like data) generated on the original 1000-tips trees with the HKY model  
__Annotations_20__ : tip annotations (protein-like data) generated on the original 1000-tips trees with the JTT model  
__scenarios_4__ : the list of true ancestral states (A, T, G or C) for all internal nodes of the original tree, which were observed in the data generation  
__scenarios_20__ : the list of true ancestral states (20 amino-acids) for all internal nodes of the original tree, which were observed in the data generation  
__PASTML-v0.6.6.1__ : PASTML source code and binary to be used in the ancestral state reconstructions  
__Seq-Gen.v.1.3.3__ : Seq-Gen source code and binary to be used fo generating *Annotations* and *scenarios* 
__Results_simulation__: Results of the accuracy of ML- and Parsimony-based methods. R commands are included to visualize the results. 
  
## Results
  
A table format of the summary statistics will be generated in each folder of *Results_simulation*
  
__Pars\_(NUC,AA)\_DOWNPASS__ : The files of 'Results.parsimony.1000.(4,20).downpass.txt' contain the statistics of the Parsimony (DOWNPASS) method. First column shows the speciation substitution rate ratios. Second column shows Brier score of the ancestral states predicted for a certain node, computed and averaged among all trees and data generated for the rate. Third columns shows the Brier score of ancestral states predicted for a certain edge. Fourth column shows the average number of states predicted for a certain node.
  
  
__NUC\_(JC,F81,HKY)__ : The files of 'Results.ml.1000.4.(JC,F81,HKY).original.txt' contain the statistics of the maximum-likelihood methods, Joint, Marginal, MAP, and MPPA, in analyses of DNA-like data. First column shows the speciation substitution rate ratios. Second-fifth columns show the node Brier scores computed for the four methods for each rate. 6th column shows the average number of states predicted by MPPA for a certain node. Remained columns show the edge Brier scores for the four methods.
  
  
__NUC\_(JC,F81,HKY)\_FLAT__ : The files of 'Results.ml.1000.4.(JC,F81,HKY).flat.txt' contain the same contents with __NUC\_(JC,F81,HKY)__ but here we show the results provided from the analyses with the flattened trees.
  
  
__AA\_(JC,F81,JTT)__ : The files of 'Results.ml.1000.4.(JC,F81,JTT).original.txt' contain the statistics of the maximum-likelihood methods, Joint, Marginal, MAP, and MPPA in analyses of Protein-like data. First column shows the speciation substitution rate ratios. Second-fifth columns show the node Brier scores computed for the four methods for each rate. 6th column shows the average number of states predicted by MPPA for a certain node. Remained columns show the edge Brier scores for the four methods.
  
  
__AA\_(JC,F81,JTT)\_FLAT__ : The files of 'Results.ml.1000.4.(JC,F81,JTT).flat.txt' contain the same contents with __AA\_(JC,F81,JTT)__ but here we show the results provided from the analyses with the flattened trees.
  
  
## To automatically analyze all data
  
```bash
./auto_analysis.sh
```
  
See our manuscript for details of the simulation procedure. A preprint will be available from bioRxiv soon.
