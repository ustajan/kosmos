# Minute-by-minute GCR proton spectral flux along a circular 500 km orbit, inclination 67 deg.
# Assumptions & simplifications (stated to be transparent):
# - Circular orbit, altitude = 500 km, inclination = 67 deg.
# - Geomagnetic latitude ~ geographic latitude (dipole approx).
# - Cutoff rigidity approximated by R_c = 14.5 * cos(latitude)^4 [GV].
# - Primary proton (nucleon) differential spectrum (unmodulated PDG-like):
#       I(E) = I0 * (E/1 GeV)^-2.7   with I0 = 1.8e4 m^-2 s^-1 sr^-1 GeV^-1
#   Converted to cm^-2 units (divide by 1e4).
# - No solar modulation, no east-west asymmetry, no penumbra — cutoff applied as a sharp threshold in rigidity.
# - Proton rigidity R [GV] from kinetic energy T [GeV]: p = sqrt(T^2 + 2 T m_p); R = p / Z (m_p=0.938 GeV/c^2, Z=1)
#
# Outputs:
# - Minute-by-minute time series for one orbit (period computed from orbital mechanics)
# - For each minute: cutoff rigidity R_c (GV), energy cutoff T_min (GeV), integrated fluxes above 0.1, 1, 10 GeV
# - Plots: R_c vs time, integrated fluxes vs time, and sample differential spectra at selected phases.
#
# Re-run the minute-by-minute GCR proton spectral flux code for a 500 km, 67 deg circular orbit.
# (Retry after earlier execution reset.)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import pi, sin, cos, asin, sqrt

# Orbit parameters
h_km = 500.0
incl_deg = 67.0
incl = np.deg2rad(incl_deg)
R_earth = 6371.0  # km
mu = 398600.4418  # km^3 / s^2
a = R_earth + h_km
period_s = 2 * pi * np.sqrt(a**3 / mu)
period_min = period_s / 60.0

print(f"Orbit period: {period_min:.2f} minutes")

# Time grid: minute-by-minute for one orbit (rounded up)
dt = 60.0  # seconds
n_steps = int(np.ceil(period_s / dt))
times_s = np.arange(n_steps) * dt  # seconds since epoch of this orbit sample
times_min = times_s / 60.0

# Energy grid for protons (GeV)
E = np.logspace(-1, 3, 300)  # 0.1 GeV to 1000 GeV
m_p = 0.938  # GeV

# PDG-like differential flux (m^-2 s^-1 sr^-1 GeV^-1) then convert to cm^-2
I0_m = 1.8e4
def pdg_spectrum(Egev):
    return I0_m * (Egev / 1.0)**-2.7  # m^-2 s^-1 sr^-1 GeV^-1

I_E_m = pdg_spectrum(E)
I_E = I_E_m / 1e4  # cm^-2 s^-1 sr^-1 GeV^-1

# Helper: rigidity R [GV] from kinetic energy T [GeV] for protons (Z=1)
def rigidity_from_T(T):
    p = np.sqrt(T**2 + 2*T*m_p)  # GeV/c
    return p  # GV (since Z=1, and p in GeV/c)

# Invert rigidity to get minimum kinetic energy for given R (numerical)
def Tmin_from_R(R):
    p = R
    return -m_p + np.sqrt(m_p**2 + p**2)

# Cutoff rigidity function (dipole approx) R_c = 14.5 cos^4(lat) GV
def Rc_from_lat(lat_rad):
    return 14.5 * (np.cos(lat_rad))**4

# Latitude as function of orbital phase: lat = arcsin(sin(i) * sin(u))
u = 2 * pi * times_s / period_s  # orbital phase angle
lat_rad = np.arcsin(np.sin(incl) * np.sin(u))

Rc = Rc_from_lat(lat_rad)  # GV
Tmin = np.array([Tmin_from_R(r) for r in Rc])  # GeV energy cutoff

# For each time, compute differential spectrum with sharp cutoff at Tmin (set E < Tmin → 0)
spectra = np.maximum(I_E * (E >= Tmin[:,None]), 0.0)  # shape (n_steps, nE)

# Integrated fluxes above thresholds
thresholds = [0.1, 1.0, 10.0]  # GeV
integrated = {th: [] for th in thresholds}
for i in range(n_steps):
    for th in thresholds:
        mask = E >= th
        val = np.trapezoid(spectra[i, mask], E[mask])  # cm^-2 s^-1 sr^-1 (since spectrum in per GeV)
        integrated[th].append(val)

df = pd.DataFrame({
    'time_min': times_min,
    'lat_deg': np.rad2deg(lat_rad),
    'Rc_GV': Rc,
    'Tmin_GeV': Tmin,
    'flux_gt_0.1_GeV': integrated[0.1],
    'flux_gt_1_GeV': integrated[1.0],
    'flux_gt_10_GeV': integrated[10.0]
})

# Display a chunk of the dataframe to the user
try:
    import caas_jupyter_tools as cjt
    cjt.display_dataframe_to_user("Minute_by_minute_GCR_flux_500km_inc67", df.head(60))
except Exception:
    # fallback: print a preview and continue (useful when running outside CAAS/Jupyter)
    print("caas_jupyter_tools not available — showing dataframe head instead:")
    n=df.__len__()
    print(f"Showing all {n} rows of the dataframe:")
    print(df.head(n).to_string())

# Plots
plt.figure(figsize=(10,4))
plt.plot(times_min, df['Rc_GV'])
plt.ylabel('Cutoff rigidity R_c (GV)')
plt.xlabel('Time (minutes)')
plt.title('Cutoff rigidity along one orbit (500 km, incl 67°)')
plt.grid(True)
plt.tight_layout()
plt.savefig('cutoff_rigidity_along_orbit.png')
plt.show()

plt.figure(figsize=(10,4))
plt.plot(times_min, df['flux_gt_0.1_GeV'], label='>0.1 GeV')
plt.plot(times_min, df['flux_gt_1_GeV'], label='>1 GeV')
plt.plot(times_min, df['flux_gt_10_GeV'], label='>10 GeV')
plt.yscale('log')
plt.xlabel('Time (minutes)')
plt.ylabel('Integrated flux (cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
plt.title('Integrated GCR proton flux along orbit (sharp cutoff applied)')
plt.legend()
plt.grid(True, which='both', ls=':')
plt.tight_layout()
plt.savefig('integrated_flux_along_orbit.png')
plt.show()

# Sample spectra at four phases: equator crossing (u=0), quarter (u=pi/4), max lat (u=pi/2), mid (u=3pi/4)
phases = [0, np.pi/4, np.pi/2, 3*np.pi/4]
phase_idx = [int(np.argmin(np.abs(u - ph))) for ph in phases]

plt.figure(figsize=(7,5))
for idx, lab in zip(phase_idx, ['Equator', 'Ascending mid', 'Max lat', 'Descending mid']):
    plt.loglog(E, spectra[idx,:], label=f'{lab} lat={df.loc[idx,"lat_deg"]:.1f}° Rc={df.loc[idx,"Rc_GV"]:.2f} GV')
plt.xlabel('Kinetic energy (GeV)')
plt.ylabel('Differential flux (cm$^{-2}$ s$^{-1}$ sr$^{-1}$ GeV$^{-1}$)')
plt.title('Sample differential GCR proton spectra at selected orbital phases')
plt.legend()
plt.grid(True, which='both', ls=':')
plt.tight_layout()
plt.savefig('sample_spectra_at_phases.png')
plt.show()

# Save csv for user
csv_path = 'gcr_flux_500km_inc67_minute_by_minute.csv'
df.to_csv(csv_path, index=False)
print(f"Saved CSV: {csv_path}")
