# generalized-method-of-moments
Scripts for analyzing data samples using the generalized method of moments.

The files contained in this folder contain the code used in the following 
publication:

Piatt, S., and Price, A. C. (2019). Analyzing dwell times with the Generalized 
Method of Moments. PloS one 14, e0197726.

The scripts and functions included here are not intended to be end-user
applications.  They are code used by us to evaluate the usefulness of the
statistal method known as the Generalized Method of Moments (GMM) as it applies
to the specific data sets discussed in our publication.  This code is intended 
for research into applications of this method and should ONLY be used by those 
who are knowledgeable BOTH in Python code AND the mathematical technique 
GMM.  We make no guarantee that this code will work for any other application.

If you wish to use this code to either reproduce or extend our results,
we recommend that you first read the publication and then contact the 
communicating author directly.

INSTRUCTIONS FOR USE

This code was written using Python 2.7.  To run everything included here, you 
will need the following packages:

numpy
scipy.optimize
cPickle
random

Download all .py and .pkl files from the GenMM project.  Rename the following
files like this (commands given in linux):

>mv gmm_NNNNNNNN.py gmm.py

>mv sim_NNNNNNNN.py sim.py

>mv cumulant_NNNNNNNN.py cumulant.py

In the above, NNNNNNNN represents a version number in the downloaded file. 
Start python and enter the following commands:

>from numpy import zeros

>import gmm

>import sim

>import cumulant as cu

>cu.initialize()

At this point you should be able to run the scripts using the following 
commands:

>run -i genData_1step_NNNNNNNN.py

>run -i genData_2step_1pass_NNNNNNNN.py

>run -i genData_2step_2pass_NNNNNNNN.py

>run -i genData_3step_NNNNNNNN.py

Again, NNNNNNNN should be replaced with a version number.  After each script 
is executed, the data will be stored in the arrays define in the scripts.  
See each script for details.
















