import os
from pastml import JOINT, MARGINAL_APPROXIMATION, MAX_POSTERIORI, JC, F81, DOWNPASS, ACCTRAN, DELTRAN, EFT

from cytopast.pastml_analyser import pastml_pipeline

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TREE_NWK = os.path.join(DATA_DIR, 'Albanian.tree.152tax.tre')
STATES_INPUT = os.path.join(DATA_DIR, 'data.txt')

if '__main__' == __name__:
    # The initial tree without ACR
    pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                    html=os.path.join(DATA_DIR, 'trees', 'Albanian_tree_initial.html'),
                    data_sep=',', verbose=True, copy_columns=['Country'], work_dir=os.path.join(DATA_DIR, 'pastml'))
    # ACR with ML methods
    for model in (F81, JC, EFT):
        for method in (JOINT, MAX_POSTERIORI, MARGINAL_APPROXIMATION):
            pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                            html_compressed=os.path.join(DATA_DIR, 'maps',
                                                         'Albanian_map_{}_{}.html'.format(method, model)),
                            html=os.path.join(DATA_DIR, 'trees', 'Albanian_tree_{}_{}.html'.format(method, model)),
                            data_sep=',', model=model, verbose=True, prediction_method=method, tip_size_threshold=1e5,
                            work_dir=os.path.join(DATA_DIR, 'pastml'))
    # ACR with parsimony
    for method in (DOWNPASS, ACCTRAN, DELTRAN):
            pastml_pipeline(data=STATES_INPUT, tree=TREE_NWK,
                            html_compressed=os.path.join(DATA_DIR, 'maps',
                                                         'Albanian_map_{}.html'.format(method)),
                            html=os.path.join(DATA_DIR, 'trees', 'Albanian_tree_{}.html'.format(method)),
                            data_sep=',', verbose=True, prediction_method=method, tip_size_threshold=1e6,
                            work_dir=os.path.join(DATA_DIR, 'pastml'))
