## Simulations of proton background in the single crystal diamond detector 

This directory contains the inputs for simulating the proton induced counts and energy deposition in a 1mm-thick diamond detector.  

Files:  

*   `detector_diamond.gdml` input to grasshopper. This models an EJ-276 surrounded on all sides by diamond.
*   `detector_flat.gdml` input to grasshopper. This models an EJ-276 covered on top and bottom by diamond.
*   `diff.AP9.output_mean_flux.AP9_i316.tsv` the spectrum of electrons.  The file name says it all.

Before running one needs to copy the `.tsv` file above to `input_spectrum.txt`.

### Running grasshopper in batch mode

In order to run a batch run, use a modification of the following bash command:  
` for (( seed=0*14; seed<60*14; seed+=14 )); do ~/git/grasshopper/exec/grass_batch.sh detector_flat.gdml diamond $seed ; done`

The script above assumes 14 cores. Replace `14` with your # of cores.  In the example above `grass_batch.sh` will run #core=14 number of simultaneous runs.  The for loop then will run 60\*14=840 runs. The way `detector_diamond.gdml` is encoded, 14 runs will reflect 20 seconds worth of incident protons.  Thus, for the example above, 60\*20=1200 seconds of exposure will be simulated. 

###Postprocessing
To postprocess the output of the simulation above, use the following bash command:  
`for ((i=0;i<840;i++)); do python ~/git/kosmos/postprocessing/sum_proton_E.py diamond_dead${i}.dat -o sum${i}.dat ; done `

The output of the above will the individual `sumX.dat` files with the total energy deposited in the internal and external (veto) detectors, with threshold limits hard-coded in `sum_proton_E.py`.  Ideally you want to see near-zero output.  Whatever survives the filtering/summing, will be the proton-induced neutron background in the detector. 

