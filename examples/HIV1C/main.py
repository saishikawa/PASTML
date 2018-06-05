import os
from pastml import *

from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE = os.path.join(DATA_DIR, 'tree.nwk')
STATES_INPUT = os.path.join(DATA_DIR, 'data.tab')

if '__main__' == __name__:
    columns = ['Location', 'RT:M184V']
    pastml_pipeline(data=STATES_INPUT, tree=TREE,
                    html=os.path.join(DATA_DIR, 'trees', 'tree_initial.html'),
                    html_compressed=os.path.join(DATA_DIR, 'maps', 'map_initial.html'),
                    verbose=True, copy_columns=columns, work_dir=os.path.join(DATA_DIR, 'pastml'),
                    name_column='Location')
    for model in (F81, JC):
        for method in (JOINT, MAX_POSTERIORI, MARGINAL_APPROXIMATION):
            pastml_pipeline(data=STATES_INPUT, tree=TREE, columns=columns,
                            html_compressed=os.path.join(DATA_DIR, 'maps', 'map_{}_{}.html'.format(model, method)),
                            html=os.path.join(DATA_DIR, 'trees', 'tree_{}_{}.html'.format(model, method)),
                            model=model, verbose=True, prediction_method=method,
                            work_dir=os.path.join(DATA_DIR, 'pastml'), name_column='Location')

    for method in (DOWNPASS, ACCTRAN, DELTRAN):
            pastml_pipeline(data=STATES_INPUT, tree=TREE, columns=columns,
                            html_compressed=os.path.join(DATA_DIR, 'maps', 'map_{}.html'.format(method)),
                            html=os.path.join(DATA_DIR, 'trees', 'tree_{}.html'.format(method)),
                            verbose=True, prediction_method=method,
                            work_dir=os.path.join(DATA_DIR, 'pastml'), name_column='Location')
