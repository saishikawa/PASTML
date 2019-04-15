# PastML C library
This repository contains C code for PastML versions till 1.0. However we have still decided to move to python3. You can find the new versions of PastML and more information at [evolbioinfo/pastml](https://github.com/evolbioinfo/pastml).

# Article

For a detailed description of PastML see Ishikawa SA, Zhukova A, Iwasaki W, Gascuel O (2019) __A Fast Likelihood Method to Reconstruct and Visualize Ancestral Scenarios__ [[bioRxiv]](https://doi.org/10.1101/379529).

# Try it online
Try it at [pastml.pasteur.fr](https://pastml.pasteur.fr)

# Input data
As an input, one needs to provide a **rooted** phylogenetical tree in [newick](https://en.wikipedia.org/wiki/Newick_format) format, and a table containing tip states.

# Run it on your computer

To run the core PastML __C library__ (__without visualisation__), follow the instructions below. For running PastML with visualisation see [evolbioinfo/pastml](https://github.com/evolbioinfo/pastml).

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

First install cmake (> ver.3.6) and GNU Scientific Library (> ver.2.3), following the instructions on [CMAKE web site](https://cmake.org/download/) and [GNU GSL web site](https://www.gnu.org/software/gsl/). GCC (> ver.4.8.5) is also required.

Then download and extract PastML, change to its folder and run:
```bash
cmake .
```
```bash
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

## Examples

See the [examples folder](https://github.com/saishikawa/PASTML/tree/master/examples) for ideas :)
