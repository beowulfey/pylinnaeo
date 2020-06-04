Quick note -- this program is still being actively developed! It's not quite feature complete, although it's close! If you like it but find something is missing or broken, please leave a note in the issues tab (or wait a few days -- I'm adding stuff pretty regularly at the moment)

Thanks! --beowulfey

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

#### Nice. What's it look like? 

Here is a screenshot as of Jun 3, 2020:

![screenshot](linnaeo/resources/images/tt-example.png "Title")


#### What's the latest? 
Here's a basic feature list so far (v0.1.6): 
* Import/export of protein sequences and alignments
* Save/Load workspaces to keep your work
* Combine sequences as needed to make new alignments using ClustalO
* Sequences and alignments can be renamed and organized in subfolders as needed
* Alignments wrap dynamically within the size of the window
* Includes a ruler and a tooltip popup on click that shows the currently highlighted residue (and the corresponding numbers for the other sequences)
* Preliminary residue theming, multiple themes available (more to come, subject to change)
* Font and point size selection for the window
* Can export the alignment window to PNG. 

### What's up next?

There is a lot I still want to do! Check out the "Issues" tab for stuff I know is broken. A quick list of features I want to add:

* Currently only supports the Clustal Omega algorithm, and calls it even for just two sequences
* Alignment order is not captured from ClustalO --> it just returns the input order, which is deceiving and confusing if you are used to the ClustalO website behavior
* I want to add better drag-and-drop functionality, such as dragging a new sequence onto an alignment window, deleting sequences from an alignment, or rerranging the order (the former 2 would be non-destructive)
* I need to make a preferences panel, with default preferences stored upon loading the app
* I want to add more color themes
* I want to make a dark theme for the app itself
* I want to incorporate a PDB viewer
* I want to show secondary structure within the alignment window.
* Residue annotations for marking comments and observations. 

#### How do I get set up? ####

## [Downloads](https://drive.google.com/drive/folders/1uk4Vd8ioxuDsuYsToDWuR-IZZBSoDJhy?usp=sharing) ##
I've attempted to build portable binaries for all three platforms. Hopefully they work, but I'm still learning this part. 

Download the correct platform, unzip, and run the linnaeo(.exe,.app) binary. Note that the executable HAS to be in the same directory as the other libraries -- it won't work otherwise! Symlinks don't work on linux (but shortcuts on Windows do work, so you can add it to your start menu). 

If you want to try building it yourself, or they don't work, here are some instructions:

###### MAC
Install homebrew if you don't have it... it will make your life significantly easier. 
Then follow these steps:
```
brew install python3
pip3 install virtualenv
brew install clustal-omega
brew install openmpi
```
Move to where you want the app to live and clone Linnaeo, and prep your env:
``` 
cd ~/devel
git clone https://github.com/beowulfey/linnaeo.git
cd linnaeo
python3 -m virtualenv venv
source venv/bin/activate
```
Install PyClustalo into it first. This is the trickiest one, and I hope it works for you. 
```
export CC=gcc-9
pip3 install --global-option=build_ext --global-option="-I/usr/local/include/clustalo" --global-option="-L/usr/local/lib" clustalo
``` 
Finally, install linnaeo, and run the program.
```
python3 setup.py build install
python3 -c 'import linnaeo; linnaeo.main()
```

###### Linux
Instructions are similar. You can try using the prebuilt wheel I have in the INSTALL folder, but they seem to have issues with compilation (the Python 3.8 version works better, so try using a venv with Py3.8). I haven't debugged this but it stopped working after I updated to Ubuntu 20.04. You can try yourself with your own compiled version of clustalo -- I'll try and update this soon with more information.

###### Windows
Install Anaconda3 and create a new environment:

```
conda create env --name linnaeo python=3.7
cd C:\Users\<You>\devel\ 	# or where ever you want it to live
git clone https://github.com/beowulfey/linnaeo.git
cd linnaeo
python setup.py build install
```
On Windows, I have a pre-compiled clustalo wheel file that seems to usually work. Please let me know if it doesn't. 


Repositories I am eternally grateful for -- they helped me get this onto windows -- and need to cite:

* [ARGTABLE2](https://github.com/jonathanmarvens/argtable2) -- for building Clustal Omega on Windows
* [Clustal Omega, adapted to use CMake (so, so grateful) from GSL Biotech](https://github.com/GSLBiotech/clustal-omega/tree/master/src)

