#!/usr/bin/env python3
"""
process_upside.py

Process 9U_upside.dat to produce filtered/grouped TSVs, compute kinematics,
and produce histogram data as NPZ and PNG.

Usage:
    python process_upside.py [input_path]

If no input_path is provided, defaults to "9U_upside.dat".
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Constants
M_N = 1.67492749804e-27  # kg
E_CHARGE = 1.602176634e-19  # J per eV
C_LIGHT = 299_792_458.0  # m/s
RNG_SEED = 42

def ensure_outdir(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)

def step_abc(input_path: Path, out_dir: Path):
    # Read whitespace-delimited input in chunks to handle large files
    COL_EVENT = "EventID"
    COL_PARTICLE = "ParticleName"
    COL_DET = "detector#"
    COL_E = "E_deposited(MeV)"
    COL_X = "x_incident"
    COL_Y = "y_incident"
    COL_Z = "z_incident"
    COL_T = "Time"

    out_a = out_dir / "a_proton_filtered_upside.tsv"
    out_b = out_dir / "b_grouped_by_event_detector_upside.tsv"
    out_c = out_dir / "c_minmax_detector_pairs_upside.tsv"
    out_d = out_dir / "d_final_metrics_upside.tsv"

    ensure_outdir(out_a)

    chunksize = 100_000
    reader = pd.read_csv(input_path, sep=r"\s+", engine="python", chunksize=chunksize)

    # Accumulator for grouped aggregations across chunks:
    # index = (EventID, detector#), columns = E_sum, x_sum, y_sum, z_sum, t_sum, cnt
    agg_acc = None
    first_chunk_written = False

    for chunk in reader:
        # (a) Filter to protons (case-insensitive)
        if COL_PARTICLE in chunk.columns:
            protons = chunk[chunk[COL_PARTICLE].str.lower() == "proton"].copy()
        else:
            protons = pd.DataFrame(columns=chunk.columns)

        # append filtered rows to out_a incrementally
        if not protons.empty:
            mode = "w" if not first_chunk_written else "a"
            header = not first_chunk_written
            protons.to_csv(out_a, sep="\t", index=False, mode=mode, header=header)
            first_chunk_written = True

        # build per-chunk group aggregates
        if not protons.empty:
            grp = (
                protons
                .groupby([COL_EVENT, COL_DET], as_index=True)
                .agg(
                    {
                        COL_E: "sum",
                        COL_X: "sum",
                        COL_Y: "sum",
                        COL_Z: "sum",
                        COL_T: "sum",
                    }
                )
            )
            # add count (number of rows per group) for mean calculation
            cnt = protons.groupby([COL_EVENT, COL_DET]).size().rename("cnt")
            grp = grp.join(cnt)

            grp.columns = ["E_sum", "x_sum", "y_sum", "z_sum", "t_sum", "cnt"]

            if agg_acc is None:
                agg_acc = grp.astype(float)
            else:
                agg_acc = agg_acc.add(grp, fill_value=0.0)

    # If no proton rows found, create empty outputs and return
    if agg_acc is None or agg_acc.shape[0] == 0:
        # create empty b,c,d files
        ensure_outdir(out_b)
        pd.DataFrame(columns=[COL_EVENT, COL_DET, "E_deposited_sum(MeV)", "x_incident_mean", "y_incident_mean", "z_incident_mean", "Time_mean"]).to_csv(out_b, sep="\t", index=False)
        ensure_outdir(out_c)
        pd.DataFrame().to_csv(out_c, sep="\t", index=False)
        ensure_outdir(out_d)
        pd.DataFrame().to_csv(out_d, sep="\t", index=False)
        return out_a, out_b, out_c, out_d, pd.DataFrame()

    # build b from aggregated sums: means = sum / cnt
    agg_acc = agg_acc.fillna(0.0)
    agg_acc["x_mean"] = agg_acc["x_sum"] / agg_acc["cnt"]
    agg_acc["y_mean"] = agg_acc["y_sum"] / agg_acc["cnt"]
    agg_acc["z_mean"] = agg_acc["z_sum"] / agg_acc["cnt"]
    agg_acc["t_mean"] = agg_acc["t_sum"] / agg_acc["cnt"]
    agg_acc = agg_acc.reset_index().rename(
        columns={
            "E_sum": "E_deposited_sum(MeV)",
            "x_mean": "x_incident_mean",
            "y_mean": "y_incident_mean",
            "z_mean": "z_incident_mean",
            "t_mean": "Time_mean",
        }
    )
    b = agg_acc[[COL_EVENT, COL_DET, "E_deposited_sum(MeV)", "x_incident_mean", "y_incident_mean", "z_incident_mean", "Time_mean"]].sort_values([COL_EVENT, COL_DET])
    ensure_outdir(out_b)
    b.to_csv(out_b, sep="\t", index=False)

    # (c) For each EventID, pick the min detector and the next (second smallest) detector
    counts = b.groupby(COL_EVENT)[COL_DET].nunique()
    eligible_events = counts[counts >= 2].index

    b_eligible = (
        b[b[COL_EVENT].isin(eligible_events)]
        .sort_values([COL_EVENT, COL_DET])
        .copy()
    )
    b_eligible["ord"] = b_eligible.groupby(COL_EVENT).cumcount()

    min_rows = b_eligible[b_eligible["ord"] == 0].copy()
    max_rows = b_eligible[b_eligible["ord"] == 1].copy()

    suffix_map = {
        "E_deposited_sum(MeV)": "E_low(MeV)",
        "x_incident_mean": "x_low",
        "y_incident_mean": "y_low",
        "z_incident_mean": "z_low",
        "Time_mean": "Time_low",
        COL_DET: "det_low",
    }
    min_rows = min_rows.rename(columns=suffix_map)

    suffix_map_hi = {
        "E_deposited_sum(MeV)": "E_high(MeV)",
        "x_incident_mean": "x_high",
        "y_incident_mean": "y_high",
        "z_incident_mean": "z_high",
        "Time_mean": "Time_high",
        COL_DET: "det_high",
    }
    max_rows = max_rows.rename(columns=suffix_map_hi)

    c = pd.merge(
        min_rows[[COL_EVENT, "det_low", "E_low(MeV)", "x_low", "y_low", "z_low", "Time_low"]],
        max_rows[[COL_EVENT, "det_high", "E_high(MeV)", "x_high", "y_high", "z_high", "Time_high"]],
        on=COL_EVENT,
        how="inner",
    )

    # Ensure det_low < det_high; if not, swap the paired columns
    swap_mask = c["det_low"] > c["det_high"]
    swap_cols = [
        ("det_low", "det_high"),
        ("E_low(MeV)", "E_high(MeV)"),
        ("x_low", "x_high"),
        ("y_low", "y_high"),
        ("z_low", "z_high"),
        ("Time_low", "Time_high"),
    ]
    for lo, hi in swap_cols:
        tmp = c.loc[swap_mask, lo].copy()
        c.loc[swap_mask, lo] = c.loc[swap_mask, hi].values
        c.loc[swap_mask, hi] = tmp.values

    ensure_outdir(out_c)
    c.to_csv(out_c, sep="\t", index=False)

    # (d) Final metrics
    d = c.copy()
    d["dTime"] = pd.to_numeric(d["Time_high"], errors="coerce") - pd.to_numeric(d["Time_low"], errors="coerce")
    for col in ["x_low", "x_high", "y_low", "y_high", "z_low", "z_high"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d["dx_incident"] = d["x_high"] - d["x_low"]
    d["dy_incident"] = d["y_high"] - d["y_low"]
    d["dz_incident"] = d["z_high"] - d["z_low"]

    d_final = d[[COL_EVENT, "det_low", "det_high", "E_low(MeV)", "E_high(MeV)", "dTime", "dx_incident", "dy_incident", "dz_incident"]].sort_values(COL_EVENT)
    ensure_outdir(out_d)
    d_final.to_csv(out_d, sep="\t", index=False)

    return out_a, out_b, out_c, out_d, d_final

def compute_kinematics(in_path: Path, out_path: Path, rng_seed: int = RNG_SEED):
    df = pd.read_csv(in_path, sep="\t")
    # Ensure numeric types
    for col in ["dx_incident", "dy_incident", "dz_incident", "dTime", "E_low(MeV)", "E_high(MeV)"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    rng = np.random.default_rng(rng_seed)
    n = len(df)

    # (a) r = sqrt(dx^2 + dy^2 + dz^2) [mm] with smearing ~10 mm (1 cm)
    dx = df["dx_incident"].to_numpy() + rng.normal(0.0, 10.0, size=n)
    dy = df["dy_incident"].to_numpy() + rng.normal(0.0, 10.0, size=n)
    dz = df["dz_incident"].to_numpy() + rng.normal(0.0, 10.0, size=n)
    df["r_mm"] = np.sqrt(dx**2 + dy**2 + dz**2)

    # (b) angle to z-axis
    rho = np.sqrt(dx**2 + dy**2)
    df["angle_to_z_rad"] = np.arctan2(rho, dz)

    # (c) velocity from dTime [ns], r [mm]
    dt_ns = df["dTime"].to_numpy() + rng.normal(0.0, 1.5, size=n)
    dt_ns_abs = np.abs(dt_ns)
    # Avoid divide-by-zero
    speed_mm_per_ns = np.where(dt_ns_abs > 0, df["r_mm"].to_numpy() / dt_ns_abs, np.nan)
    df["speed_m_per_s"] = speed_mm_per_ns * 1e6  # 1 mm/ns = 1e6 m/s

    # (d) kinetic energy (non-relativistic) in MeV
    KE_J = 0.5 * M_N * (df["speed_m_per_s"].to_numpy() ** 2)
    df["KE_MeV"] = KE_J / E_CHARGE / 1e6

    # Upside trick: confuse low/high energies and add 20% Gaussian resolution
    e_high = df["E_low(MeV)"].to_numpy() * (1.0 + rng.normal(0.0, 0.2, size=n))
    e_low = df["E_high(MeV)"].to_numpy() * (1.0 + rng.normal(0.0, 0.2, size=n))

    # Compute scattering-angle estimate: arcsin(sqrt(e/(KE+e))) clipped to [0,1]
    def safe_arcsin_arg(e, KE):
        arg = np.divide(e, (KE + e), out=np.zeros_like(e, dtype=float), where=(KE + e) > 0)
        arg = np.clip(arg, 0.0, 1.0)
        return np.sqrt(arg)

    sin_theta_low = safe_arcsin_arg(e_low, df["KE_MeV"].to_numpy())
    sin_theta_high = safe_arcsin_arg(e_high, df["KE_MeV"].to_numpy())
    theta_from_low = np.arcsin(sin_theta_low)
    theta_from_high = np.arcsin(sin_theta_high)

    # Use sign of dTime to pick which energy to use (same logic as notebook)
    dtime = df["dTime"].to_numpy()
    theta_est = np.where(dtime > 0, theta_from_low, np.where(dtime < 0, theta_from_high, np.nan))
    df["theta_error"] = theta_est - df["angle_to_z_rad"].to_numpy()

    ensure_outdir(out_path)
    df.to_csv(out_path, sep="\t", index=False)
    return df

def histogram_theta_error(df: pd.DataFrame, out_prefix: Path):
    mask = (df["E_low(MeV)"] > 0.2) & (df["E_high(MeV)"] > 0.2)
    vals = np.cos(df.loc[mask, "theta_error"].dropna().to_numpy())
    # Clip to [-1,1] to avoid nan from numerical issues
    vals = np.clip(vals, -1.0, 1.0)
    counts, edges = np.histogram(vals, bins=100, range=(-1, 1))
    total = counts.sum()
    # bins where left edge >= 0.95 and right edge <= 1.0
    bin_mask = (edges[:-1] >= 0.95) & (edges[1:] <= 1.0)
    count_bins_in_range = counts[bin_mask].sum()
    npz_path = out_prefix.with_suffix(".npz")
    np.savez(npz_path, counts=counts, edges=edges)
    # Save a PNG of the histogram
    png_path = out_prefix.with_suffix(".png")
    plt.figure(figsize=(8, 5))
    plt.hist(vals, bins=100, range=(-1, 1), log=True)
    plt.xlabel("cos(theta_error)")
    plt.ylabel("Count")
    plt.title("Histogram of cos(theta_error)")
    plt.tight_layout()
    plt.savefig(png_path, dpi=150)
    plt.close()
    return npz_path, png_path, total, count_bins_in_range, (count_bins_in_range / total if total > 0 else np.nan)

def main():
    parser = argparse.ArgumentParser(description="Process 9U_upside data and compute kinematics.")
    parser.add_argument("input", nargs="?", default="9U_upside.dat", help="Input whitespace-delimited data file (default: 9U_upside.dat)")
    parser.add_argument("--outdir", default="filtered_outputs", help="Output directory (default: filtered_outputs)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    args = parser.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_a, out_b, out_c, out_d, d_final = step_abc(input_path, out_dir)

    kin_out = out_dir / "d1_neutron_kinematics_upside.tsv"
    df_kin = compute_kinematics(out_d, kin_out)

    hist_prefix = Path("theta_error_histogram_upside")
    npz_path, png_path, total, count_in_range, frac = histogram_theta_error(df_kin, hist_prefix)

    # Minimal report prints
    if(args.verbose):
        print(f"Wrote: {out_a}")
        print(f"Wrote: {out_b}")
        print(f"Wrote: {out_c}")
        print(f"Wrote: {out_d}")
        print(f"Wrote: {kin_out}")
        print(f"Wrote: {npz_path}, {png_path}")
    print(f"Total events in histogram: {total}, counts>=0.95: {count_in_range}, fraction: {frac:.6f}")

if __name__ == "__main__":
    main()