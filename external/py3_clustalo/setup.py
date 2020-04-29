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

setup(name='clustalo',
      version='0.2.0',
      description='Python3 wrapper around libclustalo',
      author='Aaron Wolfe',
      url='https://github.com/beowolfey/clustalo-python3',
      ext_modules=[module])
