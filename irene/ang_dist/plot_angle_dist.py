'''

Plot the angular distribution of the protons from IRENE AP9 files

Date:    October 2025
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


PITCH_ANGLES = np.array([0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180])
ENERGY_GRID = np.array([
    0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1, 1.5,
    2, 3, 4, 5, 6, 7, 10, 15, 20, 30,
    40, 50, 60, 70, 100, 150, 200, 300, 400, 700,
    1200, 2000
])


def integrate_above_cut(energies, flux, threshold):
    mask = energies >= threshold
    if not np.any(mask):
        return np.zeros(flux.shape[0])

    # Compute bin edges for integration
    edges = np.zeros(len(energies) + 1)
    edges[1:-1] = 0.5 * (energies[:-1] + energies[1:])
    edges[0] = energies[0] - 0.5 * (energies[1] - energies[0])
    edges[-1] = energies[-1] + 0.5 * (energies[-1] - energies[-2])

    bin_widths = edges[1:] - edges[:-1]

    # Integrate flux × bin width (MeV)
    return np.sum(flux[:, mask] * bin_widths[mask], axis=1)


def load_ap9_directional_flux(filename, pitch_angles=PITCH_ANGLES):
    df = pd.read_csv(filename, comment='#', header=None)

    # Extract flux values
    flux_vals = df.iloc[:, 5:].to_numpy()
    n_rows, n_energy_file = flux_vals.shape

    # Energies: use the first n_energy_file from the standard grid
    energies = ENERGY_GRID[:n_energy_file]

    # Dimensions
    n_angle = len(pitch_angles)
    n_time = n_rows // n_angle

    # Construct deterministic time axis (1 min increments)
    #elapsed_minutes = np.arange(n_time)

    # Allocate array
    flux_3d = np.zeros((n_time, n_angle, n_energy_file), dtype=float)

    # Fill deterministically by row order
    for row in range(n_rows):
        ti = row // n_angle
        ai = row % n_angle
        flux_3d[ti, ai, :] = flux_vals[row, :]

    return  energies, pitch_angles, flux_3d
    #return elapsed_minutes, energies, pitch_angles, flux_3d



def plot_ap9_proton_PAD_timeavg(direction_file, lower_energy_cut=1.0):
    # Load directional flux data
    energies, pitch_angles, flux_3d = load_ap9_directional_flux(direction_file)
    #t_min, energies, pitch_angles, flux_3d = load_ap9_directional_flux(direction_file)


    n_time, n_angle, n_energy = flux_3d.shape

    # Integrate over energy for each angle and time
    J_pad_time = np.zeros((n_time, n_angle))
    for ai in range(n_angle):
        J_pad_time[:, ai] = integrate_above_cut(energies, flux_3d[:, ai, :], lower_energy_cut)

    # Time-average PAD
    J_pad_avg = np.mean(J_pad_time, axis=0)

    # Plot PAD
    plt.figure(figsize=(8, 5))
#    plt.semilogy(pitch_angles, J_pad_avg, 'o-')
    plt.plot(pitch_angles, J_pad_avg, 'o-')
    plt.yscale('linear')
    plt.xlabel('Pitch angle [deg]')
    plt.ylabel(f"Proton flux [cm$^{{-2}}$ s$^{{-1}}$ sr$^{{-1}}$]\n(E ≥ {lower_energy_cut:.2f} MeV)")
    plt.title('AP9 Proton PAD (orbit-averaged)')
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.tight_layout()

    out_png = Path(direction_file).with_suffix(f'.AP9_PAD_Ege{lower_energy_cut:.1f}MeV_timeavg.png')
    plt.savefig(out_png, dpi=300)
    plt.show()

    print(f"PAD figure saved to: {out_png}")
    print("Pitch angles:", pitch_angles)
    print("PAD flux values:", J_pad_avg)

if __name__ == "__main__":
    plot_ap9_proton_PAD_timeavg('run1.AP9.output_mean_flux.txt', lower_energy_cut=0.6)
    pass
