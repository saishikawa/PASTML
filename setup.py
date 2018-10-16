from distutils.core import setup, Extension
import subprocess
import sys

import os

if "CC" in os.environ:
    print('Using {} to compile C code.'.format(os.environ['CC']))
    os.environ["CC"] += ' -x c'

# Look for GSL
try:
    proc = subprocess.Popen(['gsl-config', '--version'], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    version = out.decode('utf-8', 'ignore').rstrip()
    print("GSL version ", version, " found.")
except:
    sys.exit("GSL not found. Please install the GNU Scientific Library (https://www.gnu.org/software/gsl).")


# the C extension module
pastml_module = Extension('pastml',
                          sources=['pastmlpymodule.c', 'runpastml.c', 'tree.c',
                                   'likelihood.c', 'marginal_likelihood.c', 'marginal_approximation.c',
                                   'states.c', 'scaling.c', 'param_minimization.c', 'logger.c', 'parsimony.c'],
                          libraries=['gsl', 'gslcblas'],
                          extra_compile_args=['-std=c11'],
                          )

setup(
    name='pastml',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='1.0.8',
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
