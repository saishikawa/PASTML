from distutils.core import setup, Extension
import subprocess
import sys

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
                          sources=['pastmlpymodule.c', 'runpastml.c', 'make_tree.c',
                                   'likelihood.c', 'marginal_likelihood.c', 'marginal_approximation.c',
                                   'output_tree.c', 'output_states.c',
                                   'scaling.c', 'param_minimization.c', 'logger.c'],
                          libraries=['gsl', 'gslcblas']
                          )

setup(
    name='pastml',
    platform=['Linux', 'Windows', 'Mac OS'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='0.5.4',
    description='Python wrapper for PASTML.',
    maintainer='Anna Zhukova',
    maintainer_email='anna.zhukova@pasteur.fr',
    url='https://github.com/saishikawa/PASTML',
    download_url='https://github.com/saishikawa/PASTML',
    keywords=['PASTML', 'phylogeny', 'ancestral state inference', 'likelihood'],
    ext_modules=[pastml_module],
    headers=['pastml.h', 'runpastml.h', 'make_tree.h',
             'likelihood.h', 'marginal_likelihood.h', 'marginal_approximation.h',
             'output_tree.h', 'output_states.h',
             'scaling.h', 'param_minimization.h', 'logger.h']
)
