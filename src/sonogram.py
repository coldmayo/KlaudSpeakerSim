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
