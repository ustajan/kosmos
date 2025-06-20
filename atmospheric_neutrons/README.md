## Simulations of the atmospheric neutrons

This directory contains the files for modeling the detector response from atmospheric neutrons.  

We use the model described in [this publication](https://doi.org/10.1007/s10686-019-09624-0), with the numerical model from [github](https://github.com/pcumani/LEOBackground/tree/master).

Below is the list of the relevant files:  

 *  `atmospheric_neutrons_-78deg_2000km_0.026percm2pers.txt` The spectrum of the neutrons, as an input to grasshopper.  The filename provides the magnetic latitude, the altitude, and the total integrated flux.  
 * `detector_diamond_neutrons.gdml` The description of the simulation. 

 