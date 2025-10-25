'''

Plot differential flux from IRENE AE9/AP9 output files and integrate spallation neutron yield
using Carpenter et al. (2012) Eq. 1 scaling.

Date:    October 2025
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def integrate_spectrum_byEq1(bin_centers, bin_contents, threshold=200):

    x = np.asarray(bin_centers)
    y = np.asarray(bin_contents)

    # sort to ensure correct integration order
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    y = y[sort_idx]

    if len(x) < 2:
        raise ValueError("Need at least two points to integrate.")

    # scale only for x >= threshold
    y_mod = np.zeros_like(y)
    mask = (x >= threshold)
    y_mod[mask] = y[mask] * 50 * (0.001 * x[mask] - 0.120)  # Carpenter scaling
    y_mod[mask] *= 50 * 20  # time (s) × area (cm²)

    # integrate using Riemann sum
    integral_value = np.sum(y_mod[1:] * np.diff(x))
    return integral_value

def detect_model_type(filename):
    with open(filename, 'r') as f:
        # Read first 5 lines just in case
        header_lines = [next(f).strip() for _ in range(5)]
    # Row index 3 (4th line)
    if len(header_lines) >= 4:
        line4 = header_lines[3].upper()
        if "AE9" in line4:
            return 'AE9'
        elif "AP9" in line4:
            return 'AP9'
    return None

def load_diff_flux(filename):
    df = pd.read_csv(filename, comment='#', header=None)
    times = df.iloc[:, 0].values
    flux = df.iloc[:, 4:].values

    # Detect model type from header
    model = detect_model_type(filename)
    if model == 'AE9':
        energy_bins = np.array([
            0.04, 0.1, 0.2, 0.3, 0.4,
            0.5, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.50, 2.75, 3.0,
            3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.5, 6.0, 6.5, 7.0, 8.5,
            10.0
        ])
    else:
        energy_bins = np.array([
            0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1, 1.5,
            2, 3, 4, 5, 6, 7, 10, 15, 20, 30,
            40, 50, 60, 70, 100, 150, 200, 300, 400, 700,
            1200, 2000
        ])

    # Match bins to number of columns in file
    n_flux_cols = flux.shape[1]
    energy_bins = energy_bins[:n_flux_cols]

    # Convert MJD to minutes since start
    elapsed_minutes = np.arange(len(times))
    return model, elapsed_minutes, energy_bins, flux

def plot_irene_flux(filename, model_type=None, n_curves=15, start_step=300, step_size=1):
    model, time, energy_bins, flux = load_diff_flux(filename)

    verbose = False
    # If user supplied model_type, override auto-detection
    if model_type is not None:
        model = model_type.upper()

    n_flux_cols = flux.shape[1]

    plt.figure(figsize=(12, 4))
    colors = plt.cm.jet(np.linspace(0, 1, n_curves))

    total_neutron_count = 0  # accumulator for average neutron rate
    for i, c in zip(range(start_step, start_step + n_curves * step_size, step_size), colors):
        if i >= len(time):
            break
        y = flux[i, :n_flux_cols]
        elapsed_minutes = i 
        spallation_neutron_count = integrate_spectrum_byEq1(
            energy_bins, y, threshold=200
        )
        total_neutron_count += spallation_neutron_count
        if i == start_step and verbose:
            print("Energy bins and flux contents at this step:")
            for idx, (e_val, f_val) in enumerate(zip(energy_bins, y)):
                print(f"  Bin {idx:2d}: E = {e_val:.4f} MeV, Flux = {f_val:.6e} cm^-2 s^-1 MeV^-1")

        plt.plot(
            energy_bins,
            y,
            color=c,
            label=f'{elapsed_minutes:.0f} min, neutrons = {spallation_neutron_count:.1e} s$^{{-1}}$'
        )

    if model == 'AE9':
        ylabel = 'electron flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]'
    else:
        ylabel = 'proton flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]'

    plt.xscale('linear')
    plt.yscale('log')
    plt.xlabel('E [MeV]')
    plt.ylabel(ylabel)
    plt.title(f'IRENE {model} model – Differential flux')
    plt.legend(fontsize=8, ncol=2)
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.tight_layout()

    out_file = Path(filename).with_suffix(f'.{model}_flux_plot.png')
    plt.savefig(out_file, dpi=300)
    plt.show()
    print(f'Plot saved to: {out_file}')

    # Print average neutron rate
    avg_neutron_rate = total_neutron_count / n_curves
    print(f'Average spallation neutron rate over plotted steps: {avg_neutron_rate:.3e} s^-1')


def integrate_above_cut(bins, flux_2d, cut):
    # bins: (n_bins,) centers in MeV
    # flux_2d: (n_times, n_bins) differential flux in cm^-2 s^-1 MeV^-1
    b = np.asarray(bins)
    y = np.asarray(flux_2d)

    # construct bin edges from centers
    edges = np.empty(b.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (b[1:] + b[:-1])
    edges[0] = b[0] - 0.5 * (b[1] - b[0])
    edges[-1] = b[-1] + 0.5 * (b[-1] - b[-2])

    # widths if no threshold
    widths = edges[1:] - edges[:-1]

    # apply threshold: effective width of each bin above 'cut'
    # For each bin i, width' = max(0, edge[i+1] - max(edge[i], cut))
    left = np.maximum(edges[:-1], cut)
    right = edges[1:]
    eff_widths = np.clip(right - left, a_min=0.0, a_max=None)

    # integrate: sum_j y[:, j] * eff_widths[j]
    return (y * eff_widths).sum(axis=1)

def plot_time_resolved_neutron_yield_from_diff(ae9_file, ap9_file, threshold=200, lower_energy_cut=1.0):

    # Load both differential flux files
    model_e, time_e, bins_e, flux_e = load_diff_flux(ae9_file)
    model_p, time_p, bins_p, flux_p = load_diff_flux(ap9_file)

    #mask_e = bins_e >= lower_energy_cut
    mask_p = bins_p >= lower_energy_cut

    neutron_yields = []
    for i in range(flux_p.shape[0]):
        neutron_count = integrate_spectrum_byEq1(
            bins_p[mask_p],
            flux_p[i, mask_p],
            threshold=threshold
        )
        neutron_yields.append(neutron_count)
    neutron_yields = np.array(neutron_yields)

    # Integrate electron and proton flux over energy bins above lower_energy_cut

    flux_e_plot = integrate_above_cut(bins_e, flux_e, lower_energy_cut)  # cm^-2 s^-1
    flux_p_plot = integrate_above_cut(bins_p, flux_p, lower_energy_cut)  # cm^-2 s^-1

#   # Alternative method using searchsorted
#    start_idx = np.searchsorted(bins_e, lower_energy_cut)
#    flux_e_plot = (np.diff(bins_e)[start_idx-1:]*flux_e[:,start_idx:]).sum(axis=1)

    fig, ax1 = plt.subplots(figsize=(10, 5))

    # Electron and proton flux time series
    ax1.plot(time_e, flux_e_plot, 'r-', 
         label=f'{model_e} e$^-$ flux≥ {lower_energy_cut:.2f} MeV')
    ax1.plot(time_p, flux_p_plot, 'b-', 
         label=f'{model_p} p$^+$ flux ≥ {lower_energy_cut:.2f} MeV')
#    ax1.set_yscale('log')
    ax1.set_xlabel('Time [min]')
    ax1.set_ylabel('Charged particle flux [cm$^{-2}$ s$^{-1}$]')
    ax1.grid(True, which='both', linestyle='--', alpha=0.6)

    # Neutron yield on secondary axis
    ax2 = ax1.twinx()
    ax2.plot(time_p, neutron_yields, 'g-', label='Spallation Neutron Yield')
 #   ax2.set_yscale('log')
    ax2.set_ylabel('Neutron yield [s$^{-1}$]')

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

    plt.tight_layout()
    out_file = Path(ae9_file).with_suffix(f'.time_resolved_neutron_cut{lower_energy_cut:.1f}MeV.png')
    plt.savefig(out_file, dpi=300)
    plt.show()

    total_neutron_count = np.trapezoid(neutron_yields, x=time_p * 60)
    print(f"Time-resolved neutron yield plot saved to {out_file}")
    print(f"Total neutron count over mission (bins ≥ {lower_energy_cut} MeV): {total_neutron_count:.3e} counts")


if __name__ == "__main__":
    plot_irene_flux('diff.AE9.output_mean_flux.txt', n_curves=15, start_step=308, step_size=1) #11,0,36
    plot_irene_flux('diff.AP9.output_mean_flux.txt', n_curves=15, start_step=308, step_size=1)

    plot_time_resolved_neutron_yield_from_diff(
    'diff.AE9.output_mean_flux.txt',
    'diff.AP9.output_mean_flux.txt',
    threshold=200,
    lower_energy_cut=1.0,
)