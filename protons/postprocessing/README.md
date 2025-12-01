This directory is for post-processing the output from the grasshopper simulations. 

* `sum_proton_energy_v3.py` This python script does the post-processing.  For every entry it calculates the light output in the detector, then groups the entries by particle type and detector #, and then sums them into a total light output.
* `QF_p.csv` proton quenching factor
* `QF_c.csv` carbon quenching factor

The quenching factor tables are read from [https://doi.org/10.1016/j.nima.2010.07.087](https://doi.org/10.1016/j.nima.2010.07.087)