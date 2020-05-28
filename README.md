# LINNAEO #

Linnaeo is a python program I made mostly as an exercise in learning how to code a GUI... but also to solve a
basic problem I had: nothing out there is very good for making, storing and viewing protein alignments.


#### Protein alignments? ####

Yes, out there in the world are countless amino acid sequences, and their sequences (and small differences) are connected to
protein functions... but it's often hard to see that. This program is intended to help with this
by meeting a few basic criteria:

1) It should be fast.
2) It should store any sequence you want. (but right now, just FASTA protein sequences, sorry)
3) It should be easy to create, store, and export alignments.

That's a basic idea. Most likely those criteria will change.

#### How do I get set up? ####

###### Linux
Create a new python3.7 environment, then compile 

###### Windows
Install Anaconda3 and create a new environment:

```
conda create env --name linnaeo python=3.7
pip install biopython, psutil, pyqt5=5.9
# If on linux:
pip install --global-option=build_ext --global-option="-I/usr/local/sci/clustalo/current/include/clustalo/" --global-option="-L/usr/local/sci/clustalo/current/lib/" clustalo
# Windows:
pip install ./INSTALL/clustalo...etc.whl

# To actually build an executable:
conda install rust, pyoxidizer
```

PyOxidizer makes use of the pyoxidizer.bzl file to build! 


I'm trying to get this into a conda package that can be installed easily (akin to PyMOL -- they figured it out somehow!). Here's hoping I am successful. If I am, I'll be sure to update that here.

Repositories I am eternal grateful for and need to cite:

* [ARGTABLE2](https://github.com/jonathanmarvens/argtable2) -- for building Clustal Omega on Windows
* [Clustal Omega, adapted to use CMake (SO, SO GRATEFUL) from GSL Biotech (You guys rock, seriously)](https://github.com/GSLBiotech/clustal-omega/tree/master/src)
*
