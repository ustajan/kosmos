# kosmos

Author: Areg Danagoulian  
License: public domain  

This repository contains the Geant4 simulation analysis for the *kosmos2553* project. The manuscript can be found at [arXiv](https://arxiv.org/abs/2512.20016).

Important directories:  

* **IRENE** Contains the IRENE AP9/AE9 model output
	* Run `IRENE/first_space_runs/plot_diff_integrate.py` to generate the plots for **Fig 1** in the manuscript.  
* **figures** Directory `figures/pixel_simple` contains the `.wrl` needed for **Fig 2(a)**.
* **protons** This directory contains the studies the proton background in the detector. Data is taken from the AP-9 model.  
	* Run `proton/postprocessing/sum_proton_E_twoveto.py` to get Fig **3(a)**.  
* **neutron_backprojection** This directory contains the studies of 9U detector package, with two detector planes working in coincidence. 
	* Run `neutron_backprojection/overlay_theta_histograms.ipynb` to get **Fig. 3(b)**.
* **time\_vs\_distance** Contains the jupyter notebook for plotting the detection time vs distance.
	* Run `time_vs_distance/time_vs_distance.ipynb` to get time vs distance plot in **Fig. 4**.  


### Miscellaneous:

* **spenvis** The older AE8/AP8 model outputs
* **electrons** This directory contains the studies of the electron backgrounds.  It simulates a detector with ~mm shielding undergoing bombardment by electrons.  Data is taken from the AE-8 model.  
* **neutrons** This directory contains the simulations for the source term in singles mode.  It takes spallation neutron data  from the Carpenter paper as an input. 
