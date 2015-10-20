## pyQbox - Input/Output-Tools for Qbox
 

PyQbox is a Python module designed for an intuitive manipulation of
[Qbox](http://qboxcode.org) input and output files. It was written with special focus on the features of [IPython](http://ipython.org) such as tab completion and easy access to help docstrings via the question mark operator. 

To use pyQbox in your Python environment you have to put this folder
in the directory of your PYTHONPATH variable. On Unix/Max, if this
variable is not set, add the following line to your .bash_profile:

> export PYTHONPATH= (path_to_your_folder)

Two IPython notebooks (and support files)  are provided as an introduction
to basic inputfile handling and outputfile parsing. Example Python scripts that use pyQbox can be found in the 'demos' directory.

To use openbabel/pybel functionality, install python-openbabel or configure your openbabel source to install python bindings.

Relies upon numpy, pybel, re, copy, ... modules


Matthew Goldey Nov 2014

Modified from pyQChem, Andreas W. Hauser and Matthew Goldey, April 2014