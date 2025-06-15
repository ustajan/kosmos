## Spenvis analysis tools

This directory contains the main tools for the analysis of the AP-8 and AE-8 outputs.  

The relevant files are the following:

 * `Spenvis_CVS_reader.ipynb` A jupyter notebook which reads in the CSV files.  It perform the following tasks:
	 * reads in the electron and proton CSV files, marked as `data_e.csv` and `data_p.csv`.  The user should update the file names, or generate links to
	 	* `data_p.csv` --> `data_spectra_orbit/spenvis_spp_361min.csv   `
	 	* `data_e.csv` --> `data_spectra_orbit/spenvis_spe_361min_minimum.csv`
	 * plots the spectra for electrons and protons
	 * generates the input files to grasshopper, which are currently marked up in the repository as 
	 	* `output_protons_316_min_51419_per_cm2_per_s.tsv` -- proton input for grasshopper simulation  , which needs to be renamed to `input_spectrum.txt` before running grasshopper using the model in `<kosmos>/protons/`
	 	* `output_electrons_316_min_54897000_per_cm2_per_s.tsv` -- electron input spectrum , which needs to be renamed to `input_spectrum.txt` before running grasshopper using the model in `<kosmos>/electrons/`
* The input files for  `Spenvis_CVS_reader.ipynb` as listed above.