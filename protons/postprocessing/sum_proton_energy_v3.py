#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/env python
"""
Sum energies of protons seen in detector 0 for each EventID,
discarding any EventID that also contains a proton in detector 1.

This is an older code which reads in the quenching factors for carbon nuclei and protons, and calculates the light output in the veto signal.  DO NOT USE.

Usage:
    python sum_proton_energy.py input_file.tsv [-o output_file.tsv]
"""

import argparse
import pandas as pd
from scipy.interpolate import interp1d

# ----------------------------------------------------------------------
def filter(input_path: str, output_path: str | None) -> None:

#    df = pd.read_csv(input_path, delim_whitespace=True)
    df = pd.read_csv(input_path, sep=r'\s+')

    # Now read in the quenching factors (QF) for protons and carbon nuclei.
    # For deuterons assume half of protons.
    # For alphas assume 1/4th of protons.
    # For all nuclei assume equal to carbon. (which is a lower limit)
    response_c = LightOutputPorC('QF_c.csv')
    response_p = LightOutputPorC('QF_p.csv')

    #y_result = interp_c.get_light(x_query)    
    # Compute weighted energy contributions
    df['weighted_E'] = 0.0

    # Protons: full energy
    mask_p = df['ParticleName'] == 'proton'
    df.loc[mask_p, 'weighted_E'] = response_p.get_light(df.loc[mask_p, 'E(MeV)'])

    # Deuterons at detector 1: half energy
    mask_d = (df['detector#'] == 1) & (df['ParticleName'] == 'deuteron')
    df.loc[mask_d, 'weighted_E'] = 0.5 * response_p.get_light(df.loc[mask_d, 'E(MeV)'])

    # Alphas at detector 1: quarter energy
    mask_a = (df['detector#'] == 1) & (df['ParticleName'].str.contains('alpha|He'))
    df.loc[mask_a, 'weighted_E'] = 0.25 * response_p.get_light(df.loc[mask_a, 'E(MeV)'])

    # anything heavier
    mask_heavy = (df['detector#'] == 1) & (df['ParticleName'].str.contains('C|B|L|N')) #for carbon, Be, B, Li, and N
    df.loc[mask_heavy, 'weighted_E'] = response_c.get_light(df.loc[mask_heavy, 'E(MeV)'])

    # Sum weighted energies by EventID and detector#
    sums = (
        df
        .groupby(['EventID', 'detector#'])['weighted_E']
        .sum()
        .reset_index()
    )

    # Pivot to have one column per detector
    pivot = (
        sums
        .pivot(index='EventID', columns='detector#', values='weighted_E')
        .fillna(0)
        .reset_index()
    )

    # Rename columns for clarity
    pivot.columns = ['EventID', 'E_det0(photons)', 'E_det1(photons)']

# **Filter** to only those events where detector 1 has E > 0.6 MeV
    #filtered = pivot[pivot['E_det1(MeV)'] > 0.6]
    filtered = pivot[
    (pivot['E_det1(photons)'] < 50) &
    (pivot['E_det0(photons)'].between(600, 100000))
    ]

    # Write to TSV

    # --- Output ----------------------------------------------------------
    if output_path:
        filtered.to_csv(output_path, sep='\t', index=False)
        print(f"Saved results to {output_path}")
    else:
        print(filtered.to_string(index=False))

# ----------------------------------------------------------------------

class LightOutputPorC:
    def __init__(self, filename):
        # Read and prepare the DataFrame
        self.df = pd.read_csv(filename) # use data from "Light output response of KamLAND liquid scintillator for
                                        #protons and 12C nuclei"
        self.df.columns = self.df.columns.str.strip()
        self.df = self.df.sort_values(by='x')
        
        # Create the interpolation function
        self.interp_func = interp1d(
            self.df['x'], 
            self.df['y'], 
            kind='linear', 
            fill_value='extrapolate'
        )

    def get_light(self, energy):
        light_output_per_mev = 1e+4 # light output for EJ-301, per MeV, for electrons
        return self.interp_func(energy)*light_output_per_mev*energy  # photons
'''
interp_c = LightOutputPorC('QF_c.csv')
interp_p = LightOutputPorC('QF_p.csv')

x_query = 7
y_result = interp_c.get_light(x_query)
print(f"y({x_query}) = {y_result}")
'''

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Sum proton energies in detector 0 per EventID, "
                    "excluding EventIDs that also have a proton in detector 1."
    )
    parser.add_argument("input_file", help="Path to the input TSV file")
    parser.add_argument(
        "-o", "--output", metavar="FILE",
        help="Write result to FILE instead of printing to stdout",
        default=None,
    )
    args = parser.parse_args()
    filter(args.input_file, args.output)

