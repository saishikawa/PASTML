# Simulation datasets and scripts
  
Materials were used in the simulation study performed by Ishikawa et al., submitted to *Mol. Biol. Evol.*
  
## Data
  
__trees_1000_original__ : the pure-birth 1000-tips trees generated for 0.1 - 10 speciation substitution rate ratios  
__trees_1000_lognormal_SD05__ : same trees with *trees_1000_original* except that their branch lengths are noised by lognormal-variates independently and randomly generated for each branch under the distribution of mean=1.0 and SD=0.5  
__Annotations_4__ : tip annotations (DNA-like data) generated on the original 1000-tips trees with the HKY model  
__Annotations_20__ : tip annotations (protein-like data) generated on the original 1000-tips trees with the JTT model  
__scenarios_4__ : the list of true ancestral states (A, T, G or C) for all internal nodes of the original tree, which were observed in the data generation  
__scenarios_20__ : the list of true ancestral states (20 amino-acids) for all internal nodes of the original tree, which were observed in the data generation  
__Seq-Gen.v.1.3.3__ : Seq-Gen source code and binary to be used fo generating *Annotations* and *scenarios*  
__Results_simulation__: Results of the accuracy of ML- and Parsimony-based methods. R commands are included in /Figures to visualize the results, which were used in our manuscript. 
  
## Results
  
A table format of summary statistics will be generated in each folder of *Results_simulation*
  
__Pars\_(NUC,AA)\_DOWNPASS__ : The files of 'Results.parsimony.1000.(4,20).downpass.txt' contain statistics of the Parsimony (DOWNPASS) method. First column shows the speciation substitution rate ratios. Second column shows Brier score of the ancestral states predicted for a certain node, computed and averaged among all trees and data generated for the rate. Third columns shows Brier score of ancestral states predicted for a certain edge. Fourth column shows the average number of states predicted for a certain node.
  
  
__NUC\_(JC,F81,HKY)__ : The files of 'Results.ml.1000.4.(JC,F81,HKY).(original,lognormal).txt' contain statistics of the maximum-likelihood methods, Joint, Marginal, MAP, and MPPA, in analyses of DNA-like data. Trees kept their original branch lengths (original) or those noised by the lognormal variates (lognormal, see our manuscript). Models of Jukes and Cantor 1969 (JC), Felsenstein 1981 (F81), and Hasegawa et al. 1985 (HKY) were applied in data analyses.
  
First column shows the speciation substitution rate ratios. Second-fifth columns show the node Brier scores computed for the four methods. 6th column shows the average number of states predicted by MPPA for a certain node. Seventh-eleventh columns contain the same information but computed for the root. Remained columns show the edge Brier scores for the four methods.
  
  
__NUC\_(JC,F81,HKY)\_LOG__ : The files of 'Results.ml.1000.4.(JC,F81,HKY).lognormal.txt' contain the same contents with __NUC\_(JC,F81,HKY)__ but here we show the results provided from the analyses with noised trees.
  
  
__AA\_(JC,F81,JTT)__ : The files of 'Results.ml.1000.4.(JC,F81,JTT).(original,lognormal).txt' contain statistics of the maximum-likelihood methods, Joint, Marginal, MAP, and MPPA in analyses of Protein-like data. Trees kept their original branch lengths (original) or those noised by the lognormal variates (lognormal, see our manuscript). JC-like and F81-like models for 20 cahracters, and Junes et al. 1992 (JTT) model were applied in data analyses.
  
First column shows the speciation substitution rate ratios. Second-fifth columns show the node Brier scores computed for the four methods. 6th column shows the average number of states predicted by MPPA for a certain node. Seventh-eleventh columns contain the same information but computed for the root. Remained columns show the edge Brier scores for the four methods.
  
  
__AA\_(JC,F81,JTT)\_LOG__ : The files of 'Results.ml.1000.4.(JC,F81,JTT).flat.txt' contain the same contents with __AA\_(JC,F81,JTT)__ but here we show the results provided from the analyses with noised trees.
  
  
## To automatically analyze and summarize data
  
```bash
./auto_analysis.sh
```
  
See our manuscript for details of the simulation procedure.
