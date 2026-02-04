import numpy as np
import matplotlib.pyplot as plt
from frb.halos.generalised_halo_profiles.generalised_halo_profile import M31_GeneralizedHaloProfile
import time

# Create M31 halo model
m31 = M31_GeneralizedHaloProfile()

# Test with 1000 impact parameters
n_points = 1000
b_impact = np.linspace(0.1, 0.99, n_points)

# Compare performance: loop vs vectorized
print(f"Testing with {n_points} impact parameters...")

# Loop version (just 100 points for comparison)
print("\nLoop version (rmax=1):")
t0 = time.time()
dms_loop = [m31.DM_from_impact_param_b(b, rmax=1).value for b in b_impact[:100]]
t_loop = time.time() - t0
print(f"  100 points: {t_loop:.2f}s ({t_loop/100*1000:.1f}ms per point)")
print(f"  Estimated for {n_points}: {t_loop/100*n_points:.1f}s")

# Vectorized version with single rmax
print("\nVectorized version (rmax=1):")
t0 = time.time()
dms_vec = m31.DM_from_impact_param_b_vectorized(b_impact, rmax=1)
t_vec = time.time() - t0
print(f"  {n_points} points: {t_vec:.2f}s ({t_vec/n_points*1000:.2f}ms per point)")
print(f"  Speedup: {(t_loop/100*n_points)/t_vec:.1f}x")

# Demonstrate with different rmax per point
print("\nVectorized version with different rmax per point:")
rmax_array = np.random.choice([1, 2, 5, 10], size=n_points)
t0 = time.time()
dms_mixed = m31.DM_from_impact_param_b_vectorized(b_impact, rmax=rmax_array)
t_mixed = time.time() - t0
print(f"  {n_points} points with mixed rmax: {t_mixed:.2f}s")

# Plot comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left plot: single rmax values
ax1 = axes[0]
for rmax in [1, 2, 5, 10]:
    dms = m31.DM_from_impact_param_b_vectorized(b_impact, rmax=rmax)
    ax1.plot(b_impact, dms.value, '-', label=f'rmax = {rmax}', alpha=0.8)
ax1.set_xlabel('b_impact (b / rvir)')
ax1.set_ylabel('DM (pc/cm³)')
ax1.set_title(f'Single rmax per curve (n={n_points})')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Right plot: mixed rmax values (color by rmax)
ax2 = axes[1]
colors = {1: 'blue', 2: 'orange', 5: 'green', 10: 'red'}
for rmax_val in [1, 2, 5, 10]:
    mask = rmax_array == rmax_val
    ax2.scatter(b_impact[mask], dms_mixed.value[mask], c=colors[rmax_val],
                label=f'rmax = {rmax_val}', alpha=0.5, s=10)
ax2.set_xlabel('b_impact (b / rvir)')
ax2.set_ylabel('DM (pc/cm³)')
ax2.set_title(f'Mixed rmax values (n={n_points})')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('dm_vs_bimpact_vectorized.png', dpi=150)
plt.show()

print(f"\nPlot saved to dm_vs_bimpact_vectorized.png")
