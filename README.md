retina-simulation
=================

Code for simulating the retina.


###LNP_model_1d
A simple LNP model that takes a 1d temporal sequence and outputs spikes.  Option to use threshold or sigmoid nonlinearity.

###LNP_model_spatiotemporal
A spatio-temporal LNP model that tiles different LNP models across the visual field in a similar way to magno- and parvocellular pathways.  Available in both MATLAB and python.  

retinaFunctions.py defines the functions that return the RGC spike times given the stimulus, RFradius, RFloc, RGCtype, and sampleRate.  Here the RGC is an LNP model and so retinaFunctions.py also defines linearKernel, threshold, etc.

visual_neuron.py defines the functions that construct the different RGCs and tile them across the stimulus.  Imports retinaFunctions.

run_*.py runs these functions with a white noise stimulus.  Imports visual_neuron.  Outputs one data file per neuron with spike times listed.



###Misc. Notes
Note that some scripts require files in the "functions" repo, so make sure to add functions to working directory.

Example:
LNP_model_1d requires "col".
