# kosmos

Author: Areg Danagoulian  
License: public domain  

This repository contains the Geant4 simulation analysis for the *kosmos2553* project. 

Important directories:  

* **protons** This directory contains the studies the proton background in the detector. Data is taken from the AP-9 model.  
* **neutron_backprojection** This directory contains the studies of 9U detector package, with two detector planes working in coincidence. 
* **IRENE** Contains the IRENE AP9/AE9 model output
* **time\_vs\_distance** Contains the jupyter notebook for plotting the detection time vs distance.

Miscellaneous:

* **spenvis** The older AE8/AP8 model outputs
* **electrons** This directory contains the studies of the electron backgrounds.  It simulates a detector with ~mm shielding undergoing bombardment by electrons.  Data is taken from the AE-8 model.  
* **neutrons** This directory contains the simulations for the source term in singles mode.  It takes spallation neutron data.   from the Carpenter paper as an input. 
