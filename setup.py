from distutils.core import setup, Extension
import subprocess
import sys

import os

os.environ["CC"] = "gcc -std=gnu99 -lgsl -lgslcblas -lm"

# Look for GSL
try:
    proc = subprocess.Popen(['gsl-config', '--version'], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    version = out.decode('utf-8').rstrip()
    GSL_V = version
    print("GSL version ", version, " found.")
except:
    sys.exit("GSL not found. Please install the GNU Scientific Library (https://www.gnu.org/software/gsl).")


# the C extension module
pastml_module = Extension('pastml',
                          sources=['pastmlpymodule.c', 'runpastml.c', 'tree.c',
                                   'likelihood.c', 'marginal_likelihood.c', 'marginal_approximation.c',
                                   'states.c', 'scaling.c', 'param_minimization.c', 'logger.c', 'parsimony.c'],
                          libraries=['gsl', 'gslcblas']
                          )

setup(
    name='pastml',
    # platform=['Linux', 'Windows', 'Mac OS'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='0.6.6',
    description='Python wrapper for PASTML.',
    maintainer='Anna Zhukova',
    maintainer_email='anna.zhukova@pasteur.fr',
    url='https://github.com/saishikawa/PASTML',
    download_url='https://github.com/saishikawa/PASTML',
    keywords=['PASTML', 'phylogeny', 'ancestral state inference', 'likelihood'],
    ext_modules=[pastml_module],
    headers=['pastml.h', 'runpastml.h', 'tree.h', 'likelihood.h', 'marginal_likelihood.h', 'marginal_approximation.h',
             'states.h', 'scaling.h', 'param_minimization.h', 'logger.h', 'parsimony.h']
)
