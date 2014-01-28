retina-simulation
=================

Code for simulating the retina.


LNP_model_1d
============
A simple LNP model that takes a 1d temporal sequence and outputs spikes.  Option to use threshold or sigmoid nonlinearity.

LNP_model_spatiotemporal
========================
A spatio-temporal LNP model that tiles different LNP models across the visual field in a similar way to magno- and parvocellular pathways.  Available in both MATLAB and python.



Misc. Notes
===========
Note that some scripts require files in the "functions" repo, so make sure to add functions to working directory.

Example:
LNP_model_1d requires "col".
