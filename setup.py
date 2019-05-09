from distutils.core import setup
from Cython.Build import cythonize

directives = {'linetrace': False, 'language_level': 3}

setup(name='bwt_based',
      ext_modules=cythonize("bwt_based.py"))
