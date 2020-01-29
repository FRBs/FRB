import numpy as np
import scipy as sp

def find_quantile(grid, distribution, quantile):
    ordering = np.argsort(grid)
    dist_interpol = sp.interpolate.interp1d(np.cumsum(distribution[ordering]), grid[ordering], bounds_error=False, fill_value='extrapolate')
    quantile_ = dist_interpol(quantile)
    return quantile_