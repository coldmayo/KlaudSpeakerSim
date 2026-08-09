import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("sonogram_data.csv", index_col=0)

angles = df.index.values.astype(float)
freqs = df.columns.values.astype(float)
spl_matrix = df.values

on_axis_spl = spl_matrix[np.argmin(np.abs(angles)), :]
spl_normalized = spl_matrix - on_axis_spl

F, A = np.meshgrid(freqs, angles)

plt.figure(figsize=(12, 6))

contour = plt.contourf(F, A, spl_normalized, levels=np.linspace(-30, 6, 37), cmap='jet', extend='both')
contour = plt.pcolormesh(F, A, spl_normalized, cmap='jet', shading='gouraud')
plt.colorbar(contour, label='Normalized SPL (dB)')

plt.xscale('log')
plt.xlim(100, 20000)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Off-Axis Angle (Degrees)')
plt.title('Loudspeaker Directivity Sonogram')
plt.grid(True, which='both', color='white', linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig('Directivity.png', dpi=300)
plt.clf()

angles_rad = np.radians(angles)
target_freqs = [500, 1000, 3000, 10000]

plt.figure(figsize=(8, 8))
ax = plt.subplot(111, projection='polar')
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

for target in target_freqs:
    idx = np.argmin(np.abs(freqs - target))

    spl_line = spl_normalized[:, idx]

    floor_db = -30
    polar_display_values = np.clip(spl_line - floor_db, 0, None)
    
    ax.plot(angles_rad, polar_display_values, label=f'{target} Hz', linewidth=2)

ax.set_rticks(np.arange(0, 36, 6))
ax.set_yticklabels([str(d + floor_db) + ' dB' for d in np.arange(0, 36, 6)])
plt.title("Loudspeaker Polar Response Patterns", pad=20)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('Polar_Response.png', dpi=300)

plt.clf()

beamwidth_6db = []

for idx in range(len(freqs)):
    spl_line = spl_normalized[:, idx]
    
    above_6db_indices = np.where(spl_line >= -3.0)[0]
    
    if len(above_6db_indices) > 0:
        min_angle = angles[above_6db_indices[0]]
        max_angle = angles[above_6db_indices[-1]]
        total_span = max_angle - min_angle
        beamwidth_6db.append(total_span)
    else:
        beamwidth_6db.append(0.0)

plt.figure(figsize=(10, 5))
plt.plot(freqs, beamwidth_6db, color='crimson', linewidth=2, label='-6dB Beamwidth')
plt.xscale('log')
plt.xlim(100, 20000)
plt.ylim(0, 180)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Total Coverage Angle (Degrees)')
plt.title('Loudspeaker -6dB Beamwidth Contour')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('Beamwidth.png', dpi=300)
