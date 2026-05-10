# Plot analytic CW solutions only (A=4 km, B=0) over 4 orbits with twin y-axes.
# Vertical black lines placed at half-periods (when radial offset = -A).
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

# Parameters (corrected chief values)
Re = 6371.0  # km
perigee_alt_chief = 2000.2  # km
apogee_alt_chief = 2006.6   # km

rp = Re + perigee_alt_chief
ra = Re + apogee_alt_chief
a_chief = 0.5 * (rp + ra)
mu = 398600.4418

# period and mean motion
T = 2.0 * math.pi * math.sqrt(a_chief**3 / mu)  # seconds
T_min = T / 60.0
n = 2.0 * math.pi / T

# analytic CW parameters
A = 4.0  # km
B = 0.0  # km

# time vector covering 4 orbits, in seconds and minutes
orbits = 4
t_sec = np.linspace(0, orbits * T, 4000)
t_min = t_sec / 60.0

# analytic solutions
x_analytic = A * np.cos(n * t_sec) + B * np.sin(n * t_sec)
y_analytic = -2.0 * A * np.sin(n * t_sec) + 2.0 * B * np.cos(n * t_sec)

# prepare plot
plt.rcParams['font.size'] = 15
fig, ax_left = plt.subplots(figsize=(11,5))

# left axis: radial x (blue)
ln1 = ax_left.plot(t_min, x_analytic, label='Radial x (analytic)', color='blue', linewidth=1.8)
ax_left.set_xlabel('Time (min)')
ax_left.set_ylabel('Radial offset x (km)', color='blue')
ax_left.tick_params(axis='y', labelcolor='blue')
ax_left.set_ylim(-1.1*A, 1.1*A)
ax_left.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.6)

# right axis: along-track y (red)
ax_right = ax_left.twinx()
ln2 = ax_right.plot(t_min, y_analytic, label='Along-track y (analytic)', color='red', linewidth=1.2)
ax_right.set_ylabel('Along-track offset y (km)', color='red')
ax_right.tick_params(axis='y', labelcolor='red')
ax_right.set_ylim(-2.2*A, 2.2*A)

# vertical black lines at every half-period where radial offset = -A
# half-periods in minutes: T_min/2, 3T_min/2, 5T_min/2, ...
half_period_lines_min = np.arange(0.5 * T_min, (orbits + 0.5) * T_min, T_min)
for tl in half_period_lines_min:
    ax_left.axvline(tl, color='black', linewidth=1.0, linestyle='-')

# Legend (combine handles)
lines = ln1 + ln2
labels = [l.get_label() for l in lines]
ax_left.legend(lines, labels, loc='upper right')

#ax_left.set_title(f'Analytic CW solution — A={A} km, B={B} km over {orbits} orbits\nChief period = {T_min:.6f} min\nBlack lines: half-periods (radial = -{A} km)')

plt.tight_layout()
png_path = "relative_motion.png"
fig.savefig(png_path, dpi=200)
plt.show()

# Save CSV timeseries for reference
#df = pd.DataFrame({"time_min": t_min, "x_km": x_analytic, "y_km": y_analytic})
#csv_path = "cw_analytic_overlay_4orbits_halfperiods_timeseries.csv"
#df.to_csv(csv_path, index=False)

print("Saved plot:", png_path)
#print("Saved timeseries CSV:", csv_path)