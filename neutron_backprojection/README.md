This directory contains grasshopper simulation input for two scenarios: 

* Two panel 9U detector assembly for neutrons incident from above: `9U.gdml`.  In order to simulate the effect of spallation in diamond, one needs to set `BeamSize` to -2 and the `BeamOffsetZ` such that the neutrons originate in the top diamond detector.
* Same, but for neutrons incident from below: `9U_upside.gdml`

Additionally, it contains the jupyter notebooks for analyzing the output.  These are:

* `filter_and_analyze.ipynb` -- for analysing the output of `9U.gdml`
* `filter_and_analyze_upside.ipynb` -- for analyzing the output of `9U_upside.gdml `
* `overlay_theta_histograms.ipynb` -- for overlaying the two plots from above.


