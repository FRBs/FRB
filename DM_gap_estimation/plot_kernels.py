import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import KernelDensity
import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

def plot_kernels():
    """Visualize the KDE kernels available in Scikit-learn"""
    fig, ax = plt.subplots(figsize=(8, 6))

    X_src = np.zeros((1, 1))
    x_grid = np.linspace(-3, 3, 1000)

    for kernel in ['gaussian', 'tophat', 'epanechnikov',
                   'exponential', 'linear', 'cosine']:
        log_dens = KernelDensity(kernel=kernel).fit(X_src).score_samples(x_grid[:, None])
        ax.plot(x_grid, np.exp(log_dens), lw=3, alpha=0.5, label=kernel)
    ax.set_ylim(0, 1.05)
    ax.set_xlim
    ax.legend(fontsize=14)
    ax.set_yticklabels([])
    ax.set_xticklabels([])

plot_kernels()
plt.show()