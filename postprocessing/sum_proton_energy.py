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
    protons = df[df['ParticleName'] == 'proton']
    sums = protons.groupby(['EventID', 'detector#'])['E(MeV)'].sum().reset_index()
    pivot = sums.pivot(index='EventID', columns='detector#', values='E(MeV)').fillna(0).reset_index()
    pivot.columns = ['EventID', 'E_det0(MeV)', 'E_det1(MeV)']

# **Filter** to only those events where detector 1 has E > 0.6 MeV
    #filtered = pivot[pivot['E_det1(MeV)'] > 0.6]
    filtered = pivot[
    (pivot['E_det1(MeV)'].astype(float) < 0.1) &
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

