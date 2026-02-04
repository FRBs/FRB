import numpy as np
import matplotlib.pyplot as plt
from frb.halos.generalised_halo_profiles.generalised_halo_profile import M31_GeneralizedHaloProfile
import time

# Create M31 halo model
m31 = M31_GeneralizedHaloProfile()

# Your use case: 100 b_impact values for each rmax from 1.0 to 10.0
b_impact = np.linspace(0.1, 0.99, 100)
rmax_values = np.arange(1.0, 10.1, 0.1)  # 1.0, 1.1, 1.2, ... 10.0

n_total = len(b_impact) * len(rmax_values)
print(f"Grid dimensions: {len(b_impact)} b_impact × {len(rmax_values)} rmax = {n_total} total calculations")

# Using DM_grid method
print("\nComputing DM grid...")
t0 = time.time()
dm_grid = m31.DM_grid(b_impact, rmax_values)
t_grid = time.time() - t0
print(f"  Time: {t_grid:.2f}s ({t_grid/n_total*1000:.2f}ms per point)")

# Compare to loop estimate
print(f"  Estimated loop time: {n_total * 0.34:.0f}s ({n_total * 0.34 / 60:.1f} min)")
print(f"  Speedup: ~{n_total * 0.34 / t_grid:.1f}x")

# Plot results
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: DM vs b_impact for selected rmax values
ax1 = axes[0]
for i, rmax in enumerate([1.0, 2.0, 5.0, 10.0]):
    idx = np.argmin(np.abs(rmax_values - rmax))
    ax1.plot(b_impact, dm_grid[idx].value, '-', label=f'rmax = {rmax}')
ax1.set_xlabel('b_impact (b / rvir)')
ax1.set_ylabel('DM (pc/cm³)')
ax1.set_title('DM vs Impact Parameter')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Right: 2D heatmap of the full grid
ax2 = axes[1]
im = ax2.imshow(dm_grid.value, aspect='auto', origin='lower',
                extent=[b_impact[0], b_impact[-1], rmax_values[0], rmax_values[-1]],
                cmap='viridis')
ax2.set_xlabel('b_impact (b / rvir)')
ax2.set_ylabel('rmax (rvir)')
ax2.set_title('DM Grid (pc/cm³)')
plt.colorbar(im, ax=ax2, label='DM (pc/cm³)')

plt.tight_layout()
plt.savefig('dm_grid.png', dpi=150)
plt.show()

print(f"\nPlot saved to dm_grid.png")
print(f"Grid shape: {dm_grid.shape}")
