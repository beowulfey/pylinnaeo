import setuptools
import sys

reqs = ['pyqt5==5.9','psutil','biopython_minimal','clustalo']
deps = []
if sys.platform == 'darwin':
    deps.append('')
elif sys.platform == 'linux':
    deps.append('.install/clustalo-0.1.2-cp37-cp37m-linux_x86_64.whl')
elif sys.platform == 'win32':
    deps.append('.install/clustalo-0.1.2-cp37-cp37m-win_amd64.whl')

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="linnaeo",
    version="0.2.0",
    author="Aaron Wolfe",
    author_email="wolfe.aarond@gmail.com",
    description="A tool for creating and categorizing protein alignments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/beowulfey/linnaeo",
    packages=setuptools.find_packages(),
    install_requires=reqs,
    dependency_links=deps,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL3",
        "Operating System :: Windows",
    ],
    python_requires='>=3.7',
)
