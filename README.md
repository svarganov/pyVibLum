pyVibLum
=========

version 1.0


Authors
========

Vsevolod D. Dergachev, Sergey A. Varganov; University of Nevada, Reno 

Liviu F. Chibotaru; Katholieke Universiteit Leuven


Description
===========

pyVibLum calculates the vibronic structure of the lanthanide emission coupled to molecular vibrations and predicts the electron-vibrational coupling strengths for specific vibrations.

User manual can be found in documentation. Calculation examples are provided.


Installation
============

pyVibLum requires Python 3. We recommend Conda installation to include NumPy and SciPy dependencies.

First, create and/or activate your 'MyEnv' conda environment:

% conda create -n MyEnv python=3.6 
% conda activate MyEnv

Second, install pyVibLum using pip. In the top directory, execute

% python -m pip install .

Now, pyVibLum is installed and accessible in the 'MyEnv' environment. To use pyVibLum, activating this environment is necessary.

To deactive conda's environment:

% conda deactivate

For more information, see documentation.

Run pyVibLum using a Python scipt (examples/start.py):

% python start.py > myoutfile.out 2>&1 


Citation
========

If you use this code, please cite:

J. Phys. Chem. Lett. 2025, 16, 2309-2313
10.1021/acs.jpclett.4c03531


Documentation
=============

UserManual1.0.pdf 


Contact
=======

svarganov@unr.edu.

