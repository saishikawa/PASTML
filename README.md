# PASTML 

__PASTML__ infers ancestral states on a phylogenetical tree with annotated tips, using maximum likelihood.
The tree with reconstructed ancestral states can then be visualised as a zoomable html map with __cytopast__.


# Run PASTML

There are 3 alternative ways to run PASTML: with [docker](https://hub.docker.com/), in python3, or in C (without visualisation).


## Run PASTML with docker
As an input, one needs to provide a phylogenetical tree in [newick](https://en.wikipedia.org/wiki/Newick_format) format,
and a table containing tip states, 
in tab-delimited (by default) or csv format (to be specified with *--data_sep ,* option).

### Basic usage
```bash
docker run -v <path_to_the_folder_containing_the_tree_and_the_annotations>:/data:rw -t evolbioinfo/pastml --tree /data/<tree_file> --data /data/<annotation_file> --columns <one_or_more_column_names> --html_compressed /data/<map_name>
```

### Example
Let's assume that the tree and annotation files are in the Downloads folder, 
and are named respectively tree.nwk and states.csv.

The states.csv is a comma-separated file, containing tip ids in the first column, 
and several named columns, including *Location*, i.e.:


Tip_id | ... | Location | ...
----- |  ----- | ----- | -----
1 | ... | Africa | ...
2 | ... | Asia | ...
3 | ... | Africa | ...
... | ... | ... | ...

To reconstruct and visualise the ancestral Location states, 
one needs to run the following command:

```bash
docker run -v ~/Downloads:/data:rw -t evolbioinfo/pastml --tree /data/tree.nwk --data /data/states.csv --data_sep , --columns Location --html_compressed /data/location_map.html
```

This will produce a file location_map.html in the Downloads folder, 
that can be viewed with a browser.


###  Help

To see advanced options, run
```bash
docker run -t evolbioinfo/pastml -h
```

### Options

```
optional arguments:
  -h, --help            show the help message and exit
  -v, --verbose         print information on the progress of the analysis

annotation-related arguments:
  -d DATA, --data DATA  the annotation file in tab/csv format with the first
                        row containing the column names.
  -s DATA_SEP, --data_sep DATA_SEP
                        the column separator for the data table. By default is
                        set to tab, i.e. for tab file. Set it to ',' if your
                        file is csv.
  -i ID_INDEX, --id_index ID_INDEX
                        the index of the column in the data table that
                        contains the tree tip names, indices start from zero
                        (by default is set to 0).
  -c [COLUMNS [COLUMNS ...]], --columns [COLUMNS [COLUMNS ...]]
                        names of the data table columns that contain states to
                        be analysed with PASTML. If neither columns nor
                        copy_columns are specified, then all columns will be
                        considered for PASTMl analysis.
  --copy_columns [COPY_COLUMNS [COPY_COLUMNS ...]]
                        names of the data table columns that contain states to
                        be copied as-is, without applying PASTML (the missing
                        states will stay unresolved).

tree-related arguments:
  -t TREE, --tree TREE  the input tree in newick format.

ancestral-state inference-related arguments:
  -m {JC,F81}, --model {JC,F81}
                        the evolutionary model to be used by PASTML, by
                        default JC.
  --prediction_method {marginal_approx,marginal,max_posteriori,joint,downpass,acctran,deltran}
                        the ancestral state prediction method to be used by
                        PASTML, by default marginal_approx.
  --work_dir WORK_DIR   the working dir for PASTML to put intermediate files
                        into (if not specified a temporary dir will be
                        created).

visualisation-related arguments:
  -n NAME_COLUMN, --name_column NAME_COLUMN
                        name of the data table column to be used for node
                        names in the compressed map visualisation(must be one
                        of those specified in columns or copy_columns if they
                        are specified).If the data table contains only one
                        column it will be used by default.
  --tip_size_threshold TIP_SIZE_THRESHOLD
                        Remove the tips of size less than the threshold-th
                        from the compressed map (set to inf to keep all tips).

output-related arguments:
  -o OUT_DATA, --out_data OUT_DATA
                        the output annotation file with the states inferred by
                        PASTML.
  -p HTML_COMPRESSED, --html_compressed HTML_COMPRESSED
                        the output summary map visualisation file (html).
  -l HTML, --html HTML  the output tree visualisation file (html).
```

## Run PASTML in python3

We strongly recommend installing pastml for python via [conda](https://conda.io/docs/). 

Once you have conda installed create an environment for pastml with python3, gcc and gsl:

```bash
conda create --name pastml python=3 gcc gsl
```

Then activate it:
```bash
source activate pastml
```

Then install cytopast (a python library that adds visualisation to PASTML) in it:

```bash
pip install cytopast
```

For Windows users we recommend to use Cygwin environment, installing gcc (> 7.3.0), gsl (> 2.3), python and pip (> 3.6).
Then install pastml and cytopast:

```bash
pip3.6 install pastml
pip3.6 install cytopast
```

### Basic usage in a command line
If you installed cytopast via conda, do not forget to activate the dedicated environment, e.g.

```bash
source activate cytopast
```

To run cytopast:

```bash
cytopast --tree <path/to/tree_file.nwk> --data <path/to/annotation_file.tab> --columns <one_or_more_column_names> --html_compressed <path/to/output/map.html>
```

To see advanced options, run:
```bash
cytopast -h
```

### Basic usage in python3
```python
from cytopast.pastml_analyser import pastml_pipeline

# Path to the table containing tip/node annotations, in csv or tab format
data = "/path/to/the/table/eg/data.csv"

# Path to the tree in newick format
tree = "/path/to/the/tree/eg/tree.nwk"

# Columns present in the annotation table,
# for which we want to reconstruct ancestral states
columns = ['Location', 'Resistant_or_not']

# Columns present in the annotation table,
# for which we want to copy existing annotations from the annotation table,
# without inferring ancestral states
copy_columns = ['Sex']

# Path to the output compressed map visualisation
html_compressed = "/path/to/the/future/map/eg/map.html"

# Path to the output tree visualisation
html = "/path/to/the/future/tree/visualisation/eg/tree.html"

pastml_pipeline(data=data, data_sep=',', columns=columns, name_column='Location',
                tree=tree,
                html_compressed=html_compressed, html=html, 
                verbose=True)
```

## Run PASTML in C (without visualisation)

### Installation

First install [GNU GSL](https://www.gnu.org/software/gsl/).
Then run:
```bash
cmake
make
```


### Basic usage in a command line
```bash
pastml -t <path/to/tree_file.nwk> -a <path/to/annotation_file.csv> 
```

To see advanced options, run
```bash
cytopast -h
```

### Options
```
required arguments:
   -a ANNOTATION_FILE                  path to the annotation file containing tip states (in csv format: tip_id,state.)
   -t TREE_NWK                         path to the tree file (in newick format)

optional arguments:
   -o OUTPUT_ANNOTATION_CSV            path where the output annotation file containing node states will be created (in csv format)
   -n OUTPUT_TREE_NWK                  path where the output tree file will be created (in newick format)
   -r OUTPUT_PARAMETERS_CSV            path where the output parameters file will be created (in csv format)
   -m MODEL                            state evolution model for max likelihood prediction methods: "JC" (default) or "F81"
   -p PREDICTION_METHOD                ancestral state prediction method: "marginal_approx" (default), "marginal", "max_posteriori", "joint", "downpass", "acctran", or "deltran"
("marginal_approx", "marginal", "max_posteriori", and "joint" are max likelihood methods, while "downpass", "acctran", and "deltran" are parsimonious ones)
   -q                                  quiet, do not print progress information
```
