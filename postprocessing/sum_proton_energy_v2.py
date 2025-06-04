#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/env python
"""
Sum energies of protons seen in detector 0 for each EventID,
discarding any EventID that also contains a proton in detector 1.

Usage:
    python sum_proton_energy.py input_file.tsv [-o output_file.tsv]
"""

import argparse
import pandas as pd

# ----------------------------------------------------------------------
def filter(input_path: str, output_path: str | None) -> None:

#    df = pd.read_csv(input_path, delim_whitespace=True)
    df = pd.read_csv(input_path, sep=r'\s+')

    # Compute weighted energy contributions
    df['weighted_E'] = 0.0

    # Protons: full energy
    mask_p = df['ParticleName'] == 'proton'
    df.loc[mask_p, 'weighted_E'] = df.loc[mask_p, 'E(MeV)']

    # Deuterons at detector 1: half energy
    mask_d = (df['detector#'] == 1) & (df['ParticleName'] == 'deuteron')
    df.loc[mask_d, 'weighted_E'] = 0.5 * df.loc[mask_d, 'E(MeV)']

    # Alphas at detector 1: quarter energy
    mask_a = (df['detector#'] == 1) & (df['ParticleName'] == 'alpha')
    df.loc[mask_a, 'weighted_E'] = 0.25 * df.loc[mask_a, 'E(MeV)']

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
    pivot.columns = ['EventID', 'E_det0(MeV)', 'E_det1(MeV)']

# **Filter** to only those events where detector 1 has E > 0.6 MeV
    #filtered = pivot[pivot['E_det1(MeV)'] > 0.6]
    filtered = pivot[
    (pivot['E_det1(MeV)'] < 0.1) &
    (pivot['E_det0(MeV)'].between(0.6, 10))
    ]

    # Write to TSV

    # --- Output ----------------------------------------------------------
    if output_path:
        filtered.to_csv(output_path, sep='\t', index=False)
        print(f"Saved results to {output_path}")
    else:
        print(filtered.to_string(index=False))

# ----------------------------------------------------------------------


# In[7]:


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

