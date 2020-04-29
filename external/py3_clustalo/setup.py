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
                   include_dirs=['/usr/local/sci/clustalo/current/include/clustalo', '/usr/include/clustalo', '/usr/local/include/clustalo'],
                   libraries=libraries,
                   library_dirs=['/usr/local/sci/clustalo/current/lib'],
                   extra_compile_args=extra_compile_args)

setup(name='clustalo',
      version='0.1.2',
      description='Python wrapper around libclustalo',
      author='Joshua Ma',
      author_email='josh@benchling.com',
      url='https://github.com/benchling/clustalo-python',
      ext_modules=[module])
