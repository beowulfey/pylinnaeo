from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.
buildOptions = dict(packages = [], excludes = [])

import sys
base = 'Win32GUI' if sys.platform=='win32' else None

executables = [
    Executable('.\\linnaeo\\__init__.py', base=base, targetName = 'linnaeo')
]

setup(name='linnaeo',
      version = '0.2.0',
      description = 'Protein Alignments',
      options = dict(build_exe = buildOptions),
      executables = executables)
