import setuptools
from distutils.core import setup, Extension

import os

OPENMP_DISABLED = os.environ.get('OPENMP_DISABLED', False)
libraries=['clustalo', 'stdc++']
extra_compile_args = []

if not OPENMP_DISABLED:
    libraries.append('gomp')
    extra_compile_args.append('-fopenmp')

module = Extension('clustalo',
                   sources = ['clustalo.c'],
                   include_dirs=['./clustalo/include/clustalo'],
                   libraries=libraries,
                   library_dirs=['./clustalo/lib'],
                   extra_compile_args=extra_compile_args)

setup(name='py3-clustalo',
      version='0.2.0',
      description='Python3 wrapper around libclustalo',
      author='Aaron Wolfe',
      email='wolfe.aarond@gmail.com',
      url='https://github.com/beowulfey/clustalo-python3',
      packages=setuptools.find_packages(),
      ext_modules=[module],
      python_requires='>=3.8'
      )
