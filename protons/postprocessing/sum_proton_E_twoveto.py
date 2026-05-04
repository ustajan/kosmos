#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/env python
"""
Sum energies of protons seen in detector 0 for each EventID,
discarding any EventID that also contains a proton in detector 1.

Usage:
    python sum_proton_energy.py input_file.tsv [-o output_file.tsv]

This code is a modification of sum_proton_E.py to implement a two-veto system:
- Events are vetoed if detector 1 (downstream) records a proton energy > 0.2 MeV
  OR detector 2 (upstream) records a proton energy > 0.6 MeV.
This logic is based on ../proton/detector_flat.gdml geometries.

"""

import argparse
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ----------------------------------------------------------------------
def filter(input_path: str, output_path: str | None) -> None:

#    df = pd.read_csv(input_path, delim_whitespace=True)
    df = pd.read_csv(input_path, sep=r'\s+')

    # Sum weighted energies by EventID and detector#
    sums = (
        df
        .groupby(['EventID', 'detector#'])['E(MeV)']
        .sum()
        .reset_index()
    )

    # Pivot to have one column per detector
    pivot = (
        sums
        .pivot(index='EventID', columns='detector#', values='E(MeV)')
        .fillna(0)
        .reset_index()
    )

    # Rename columns for clarity
    pivot.columns = ['EventID', 'E_det0', 'E_det1', 'E_det2']

# **Filter** to only those events where detector 1 has E > 0.6 MeV
    #filtered = pivot[pivot['E_det1(MeV)'] > 0.6]
    filtered_veto = pivot[~(
    ((pivot['E_det1'] > 0.2) | (pivot['E_det2'] > 0.6) ) 
     #this is the energy deposition below which we have very few electrons
#    |
#    ((pivot['E_det2'] > 2) & (pivot['E_det1'] <0.2) ) #this condition comes from the assumption of a range of 1.352cm for the proton not to reach det1. E<35MeV, dE/dx>14.74, deltaE>5.18 MeV
    )]
    filtered = filtered_veto[filtered_veto['E_det0'].between(0.2, 10)]
    number_of_vetos = pivot[((pivot['E_det1'] > 0.2) | (pivot['E_det2'] > 0.6) )].__len__()
#    number_of_vetos += pivot[((pivot['E_det2'] > 4) & (pivot['E_det1'] <0.2) )].__len__()
    print(f"Number of vetoed EventIDs (E_det1 > 0.2 MeV or E_det2 > 0.6 MeV) / total events: {number_of_vetos} / {pivot.__len__()}")
    # Write to TSV

    # --- Output ----------------------------------------------------------
    if output_path:
        filtered.to_csv(output_path, sep='\t', index=False)
        print(f"Saved results to {output_path}")
    else:
        print(filtered.to_string(index=False))

    # Python
    import matplotlib.pyplot as plt
    import numpy as np

    # compute sum and ensure numeric
    sum12 = (pivot['E_det1'].astype(float) + pivot['E_det2'].astype(float)).dropna()
    x = sum12.values
    y = pivot.loc[sum12.index, 'E_det0'].astype(float).values

    plt.figure(figsize=(7,6))
    # 2D histogram (density) of veto energy vs E_det0
    from matplotlib.colors import LogNorm
    # use fixed ranges so the histogram aligns with the rectangle and axis limits
    H, xedges, yedges, im = plt.hist2d(
        x, y,
        bins=[200, 200],
        range=[[0, 6], [0, 20]],
        cmap='viridis',
        norm=LogNorm(vmin=1)
    )
    # get current axes so colorbar can be anchored to the plot (flush)
    ax = plt.gca()
    # compute bin sizes and add colorbar with label including MeV units
    bin_x = xedges[1] - xedges[0]
    bin_y = yedges[1] - yedges[0]
    #plt.colorbar(im, label=f'Counts/{bin_x:.6g} MeV/{bin_y:.6g} MeV')
    # use a small pad so the colorbar sits close/flush to the axes
    plt.colorbar(im, ax=ax, label=f'Counts/{bin_x:.6g} MeV/{bin_y:.6g} MeV', pad=0.01, fraction=0.046)
    plt.xlabel('Veto deposited energy [MeV]',fontsize=14)
    plt.ylabel(r'$E_\mathrm{EJ276}$ [MeV]',fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlim(0, 6)
    plt.ylim(0, 20)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()

    # get current axes (or use the ax returned by subplots)
    ax = plt.gca()

    # rectangle at x=[0,0.6], y=[0.5,10], blue and semi-transparent
    rect = patches.Rectangle((0.0, 0.2), 0.6, 10.0 - 0.2,
                             facecolor='blue', alpha=0.3, edgecolor='none', zorder=0)
    ax.add_patch(rect)

    # ensure scatter is drawn above the rectangle (if needed)
    # ax.scatter(x, y, s=12, alpha=0.6, edgecolors='none', zorder=2)
    plt.savefig('E_det0_vs_E_veto_twoveto.pdf')
    plt.show()
    # save raw values and histogram arrays for other notebooks
#    np.savez('E_det1_plus_E_det2_hist.npz', values=vals, counts=counts, edges=edges)

    # prepare aligned numeric arrays
    df_xy = pivot[['E_det1', 'E_det2']].astype(float).dropna()
    x = df_xy['E_det1'].values   # x-axis
    y = df_xy['E_det2'].values   # y-axis

    # 2D histogram
    plt.figure(figsize=(7,6))
    H, xedges, yedges, img = plt.hist2d(
        x, y,
        bins=100,
        range=[[0, np.percentile(x, 99.9)], [0, np.percentile(y, 99.9)]],  # adjust as needed
        cmap='viridis',
        norm=LogNorm(vmin=1)  # remove or change for linear scale
    )
    plt.colorbar(img, label='Counts')
    plt.xlabel('E_det1 [MeV]')
    plt.ylabel('E_det2 [MeV]')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig('E_det2_vs_E_det1_hist2d.pdf')
    plt.show()
# ----------------------------------------------------------------------



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

