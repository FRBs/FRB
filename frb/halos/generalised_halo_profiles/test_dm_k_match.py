import numpy as np
import matplotlib.pyplot as plt
from frb.halos.generalised_halo_profiles.generalised_halo_profile import M31_GeneralizedHaloProfile

# Impact parameters
b_impact = np.linspace(0.1, 0.99, 50)

# k values (which also serve as rmax)
k_values = [1, 2, 5, 10]

plt.figure(figsize=(10, 6))

for k in k_values:
    # Create model with k parameter matching rmax
    m31 = M31_GeneralizedHaloProfile(log_Mhalo=12.18, k=k)

    dms = []
    for b in b_impact:
        # Use rmax = k (same as the model's k parameter)
        dm = m31.DM_from_impact_param_b(b, rmax=k)
        dms.append(dm.value)

    plt.plot(b_impact, dms, 'o-', label=f'k = rmax = {k}', markersize=4)

plt.xlabel('b_impact (b / rvir)')
plt.ylabel('DM (pc/cm³)')
plt.title('DM vs Impact Parameter (k parameter matched to rmax)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('dm_vs_bimpact_k_matched.png', dpi=150)
plt.show()

print("Plot saved to dm_vs_bimpact_k_matched.png")
