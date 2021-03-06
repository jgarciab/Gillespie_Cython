import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#Call: python setup.py build_ext --inplace

ext_modules = [Extension("loopGillespie", ["loopGillespie.pyx"], include_dirs=[numpy.get_include()])]

setup(
    name = 'loopGillespie',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
