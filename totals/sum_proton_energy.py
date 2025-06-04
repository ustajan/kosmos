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
    """
    Parameters
    ----------
    input_path : str
        Path to the input tab-separated file.
    output_path : str | None
        Where to write the result (tab-separated).  If None, prints to stdout.
    """

    # Your file has no header row, so give the columns temporary names.
    # 1-based positions ➜ 0-based indices (2nd col → index 1, etc.).
    columns = [
        "col0",           # 1st column (unused here)
        "Energy",         # 2nd column
        "EventID",        # 3rd column
        "ParticleName",   # 4th column
        "col4", "col5",   # 5th–6th columns (unused)
        "DetectorNumber", # 7th column
    ]

    df = pd.read_csv(input_path, sep="\t", header=1, names=columns) #replace header=None if there is no header

    # --- Identify the rows we care about ---------------------------------
    in_det0 = (df["ParticleName"] == "proton") & (df["DetectorNumber"] == 0) & (df["Energy"] >= 0.6) & (df["Energy"] <= 10)
    in_det1 = (df["ParticleName"] == "proton") & (df["DetectorNumber"] == 1) & (df["Energy"] >= 0.3)

    # EventIDs that contain a proton in detector 1 – we will drop these.
    exclude_ids = set(df.loc[in_det1, "EventID"])

    # Keep only proton/detector 0 rows whose EventID is not in exclude_ids.
    df_det0 = df.loc[in_det0 & ~df["EventID"].isin(exclude_ids)]
    df_det1 = df.loc[in_det1 & ~df["EventID"].isin(exclude_ids)]
    

    # --- Sum the energies per EventID ------------------------------------
    summed = (
        df_det0
        .groupby("EventID", as_index=False)["Energy"]
        .sum()
        .rename(columns={"Energy": "TotalEnergy_0"})
    )

    # --- Output ----------------------------------------------------------
    if output_path:
        summed.to_csv(output_path, sep="\t", index=False)
    else:
        print(summed.to_string(index=False))


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

