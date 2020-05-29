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

I've attempted to build portable binaries for all three platforms. [They can be found here.](https://drive.google.com/drive/folders/1uk4Vd8ioxuDsuYsToDWuR-IZZBSoDJhy?usp=sharing)

Download the correct platform, unzip, and run the linnaeo(.exe,.app) binary. Note that the executable HAS to be in the same directory as the other libraries -- it won't work otherwise! Symlinks don't work on linux (but shortcuts on Windows do work, so you can add it to your start menu). 

NOTE: v0.1.5 (the current version) has a debug module in place that preloads with 2 sequences already. You can delete those, but I forgot to, so sorry about that. I promise the next version won't have that. 

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
Instructions are similar. You can try using the prebuilt wheel I have in the INSTALL folder, but they seem to have issues with compilation. I haven't debugged this but you can try yourself with your own compiled version of clustalo. I'll try and update this soon.

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
