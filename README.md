# PastML 
__PastML__ infers ancestral characters on a rooted phylogenetic tree with annotated tips, using maximum likelihood or parsimony.
The tree with reconstructed ancestral states can then be visualised as a zoomable html map with [__cytopast__](https://github.com/evolbioinfo/cytopast).

__PastML__ is a core C library that provides fast Ancestor Character Reconstruction (ACR) methods without visualisation.
For using __PastML with visualisation__ please refer to [__cytopast__ web page](https://github.com/evolbioinfo/cytopast).

# Input data
As an input, one needs to provide a **rooted** phylogenetical tree in [newick](https://en.wikipedia.org/wiki/Newick_format) format,
and a table containing tip states.

# Try it online
Try it at [pastml.pasteur.fr](https://pastml.pasteur.fr)


# Run it on your computer

If you want to run PastML __with visualisation__, please refer to the [__cytopast__ web page](https://github.com/evolbioinfo/cytopast). 
Two ways are described there: via __[docker](https://hub.docker.com/)__, and in __python3/command line__.

If you want to run the core PastML __C library__ (__without visualisation__), follow the instructions below.


## Input data example
Let's assume that we have a __rooted__ tree with n tips (in [newick](https://en.wikipedia.org/wiki/Newick_format) format), 
and a tip annotation file, and that they are in the Downloads folder, 
named respectively _tree.nwk_ and _states.csv_.

The _states.csv_ is a __comma-separated file__, containing tip ids (the same as in the tree) in the first column, 
and tip annotations in the second column. 
It must __not__ contain a header line, and it must contain an annotation for __each__ tip in the tree 
(possibly empty if the annotation is unknown), e.g.:

```text
tip_1,Africa
tip_2,Asia
tip_3,
...
tip_n,Africa
```

Note, that here we do not know the annotation for tip_3, so we left it blank.

## PastML C library Installation

First install GNU Scientific Library, following the instructions on [GNU GSL web site](https://www.gnu.org/software/gsl/).

Then download and extract PastML, change to its folder and run:
```bash
cmake
make
```

## Basic usage in a command line
```bash
pastml -t <path/to/tree_file.nwk> -a <path/to/annotation_file.csv> 
```

For the data example above, the command becomes:
```bash
pastml -t ~/Downloads/tree.nwk -a ~/Downloads/states.csv
```

## Help

To see advanced options (including the choice of evolutionary model and ancestral character reconstruction method), run:
```bash
pastml -h
```
