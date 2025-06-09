## Simulations of electron signal in the single crystal diamond detector 

This directory contains the inputs for simulating the electron counts and energy deposition in a 1mm-thick diamond detector.  

Files:  

*   `beta_diamond.gdml` input to grasshopper.  
*   `output_electrons_316_min_54897000_per_cm2_per_s.tsv` the spectrum of electrons.  The file name says it all.
*   `analysis.sh` -- simple shell script which determines the counts in the detector, the energy deposition, as well as the period / waiting time between multiple hits.  If you run it without any arguments it will print out the help with instructions.

Before running one needs to copy the `.tsv` file above to `input_spectrum.txt`.