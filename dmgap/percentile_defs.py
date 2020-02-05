import numpy as np
import scipy as sp

def find_percentile(grid, distribution, percentile):
    ordering = np.argsort(grid)
    dist_interpol = sp.interpolate.interp1d(np.cumsum(distribution[ordering]), grid[ordering], bounds_error=False, fill_value='extrapolate')
    percentile_ = dist_interpol(percentile)
    return percentile_