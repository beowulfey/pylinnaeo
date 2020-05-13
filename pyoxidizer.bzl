def make_dist():
    #return default_python_distribution()
    return default_python_distribution(flavor="standalone_dynamic")

def make_exe(dist):
    python_config = PythonInterpreterConfig(
    #     filesystem_importer=True,
    #     verbose=1,
          run_eval="import linnaeo; linnaeo.main()",
    #      run_module="linnaeo",
          sys_paths=["$ORIGIN/lib"],
    #     run_noop=False,
    #     run_repl=True,
    )

    exe = dist.to_python_executable(
        name="linnaeo",
        config=python_config,
        extension_module_filter='all',
        include_sources=True,
        include_test=False,
        #resources_policy="prefer-in-memory-fallback-filesystem-relative:lib"
    )

    for resource in dist.pip_install([
      'C:\\Users\\wolfe\\devel\\linnaeo\\install\\clustalo-0.1.2-cp37-cp37m-win_amd64.whl',
      'C:\\Users\\wolfe\\devel\\linnaeo\\install\\biopython_minimal-1.77.dev0-py3-none-any.whl',
      'C:\\Users\\wolfe\\devel\\linnaeo\\dist\\linnaeo-0.2.0-py3-none-any.whl',
      ]):
    #  'C:\\Users\\wolfe\\devel\\linnaeo\\dist\\linnaeo-0.2.0-py3-none-any.whl']):
        exe.add_in_memory_python_resource(resource)
    return exe

def make_embedded_resources(exe):
    return exe.to_embedded_resources()

def make_install(dist, exe):
    files = FileManifest()
    files.add_python_resource(".", exe)
    files.add_python_resources("lib",dist.pip_install(['pyqt5==5.9', 'psutil',
    #'C:\\Users\\wolfe\\devel\\biopy-minimal\\dist\\biopython_minimal-1.77.dev0-py3-none-any.whl',
    #'C:\\Users\\wolfe\\devel\\linnaeo\\install\\clustalo-0.1.2-cp37-cp37m-win_amd64.whl',
    #'C:\\Users\\wolfe\\devel\\linnaeo\\dist\\linnaeo-0.2.0-py3-none-any.whl'
    ]))
    return files

register_target("dist", make_dist)
register_target("exe", make_exe, depends=["dist"], default=True)
register_target("resources", make_embedded_resources, depends=["exe"], default_build_script=True)
register_target("install", make_install, depends=["dist","exe"], default=True)

resolve_targets()

# END OF COMMON USER-ADJUSTED SETTINGS.
#
# Everything below this is typically managed by PyOxidizer and doesn't need
# to be updated by people.

PYOXIDIZER_VERSION = "0.7.0"
PYOXIDIZER_COMMIT = "UNKNOWN"
