## Neutron detector's efficiency determination

This directory contains the files necessary for running simulations to determine the efficiency of the neutron detector. 

Relevant files:  

* `detector_diamond_neutrons.gdml` The input to grasshopper.  It assumes a source at the distance of 50cm from the detector. 
* `neutron_spallation.tsv` The neutron spallation spectrum, to be read into grasshopper. As always,this file needs to be copied to `input_spectrum.txt`.

Having ran 14 simulations with different seeds, with outputs named `detectorX.dat`, output of the simulation then needs to be post-processed using the usual script, as follows:  
`for ((i=0;i<14;i++)); do python ~/git/kosmos/postprocessing/sum_proton_E.py diamond${i}.dat -o sum${i}.dat ; done`

The final neutron count is determined by simply counting the lines in the output files (minus headers). 