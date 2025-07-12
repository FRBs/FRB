from importlib import resources
import os

import numpy as np
import matplotlib.pyplot as plt

import healpy as hp

from functools import partial
from concurrent.futures import ProcessPoolExecutor

from tqdm import tqdm

from ne2001 import density

from IPython import embed


def get_current_mapfile():
    print("Using the default NE2001 DM HEALPix map file.")
    """
    Get the path to the default NE2001 DM HEALPix map file.
    
    Returns
    -------
    str
        The path to the NE2001 DM HEALPix map file.
    """
    return os.path.join(
        resources.files('frb'), 'dm_ism_models','ne2001_dm_healpix_map.fits'
    )


def calc_ismDM_item(item,ne=None):
    l, b = item
    ismDM = ne.DM(l, b, 100.)
    return ismDM.value

def create_ne2001_dm_healpix_map(nside=64, n_cores=15):
    """
    Create a HEALPix map of dispersion measure (DM) values using the NE2001 model.

    Parameters
    ----------
    nside : int
        The HEALPix resolution parameter. The number of pixels will be 12 * nside^2.

    Returns
    -------
    None
    """
    # Create a HEALPix map of DM values using the NE2001 model
    npix = hp.nside2npix(nside)
    print(f"Number of pixels: {npix}")
    dm_map = np.zeros(npix)


    # Create the items
    items = []
    for pix in range(npix):
        theta, phi = hp.pix2ang(nside, pix)  # HEALPix coordinates
        l = np.degrees(phi)  # Galactic longitude
        b = 90 - np.degrees(theta)  # Galactic latitude
        items.append([l,b])

    # Run
    ne = density.ElectronDensity()
    map_fn = partial(calc_ismDM_item, ne=ne)

    # Multi-process time
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        chunksize = len(items) // n_cores if len(items) // n_cores > 0 else 1
        answers = list(tqdm(executor.map(map_fn, items,
                                            chunksize=chunksize), 
                            total=len(items)))
    dm_map = np.array(answers)

    # Save the map to a FITS file for future use
    save_path = get_current_mapfile()
    hp.write_map(save_path, dm_map, overwrite=True)

    print(f"HEALPix map created and saved as {save_path}")


def get_dm_map(mapfile=None):
    """
    Read a HEALPix map of dispersion measure (DM) values created using the NE2001 model.

    Parameters
    ----------
    mapfile : str
        The path to the FITS file containing the HEALPix map of DM values.

    Returns
    -------
    dm_map : array
        The DM map read from the FITS file.
    """
    if mapfile is None:
        mapfile = get_current_mapfile()
    # Read the HEALPix map of DM values
    dm_map = hp.read_map(mapfile)

    return dm_map

def grab_dm_ism_from_healpix_map(l,b, dm_map):
    """
    Get the dispersion measure (DM) value from a HEALPix map of DM values.

    Parameters
    ----------
    l : float
        Galactic longitude in degrees.
    b : float
        Galactic latitude in degrees.
    dm_map : array
        The HEALPix map of DM values

    Returns
    -------
    dm_ism : float
        The DM value at the given coordinates.
    """

    # Get the pixel corresponding to the given coordinates
    theta = np.radians(90 - b)  # Galactic latitude
    phi = np.radians(l)  # Galactic longitude
    pix = hp.ang2pix(nside=64, theta=theta, phi=phi)

    # Get the DM value at the given coordinates
    dm_ism = dm_map[pix]

    return dm_ism

def plot_mollwiede_view_dm_ism(mapfile=None,
                               title=r'$DM_{ISM}$ Map',
                               min=0, max=1000):
    """
    Plot a Mollweide view of the HEALPix map of dispersion measure (DM) values.

    Parameters
    ----------
    mapfile : str
        The path to the FITS file containing the HEALPix map of DM values.

    Returns
    -------
    None
    """
    if mapfile is None:
        mapfile = get_current_mapfile()
        
    # Load the HEALPix map of DM values
    dm_map = hp.read_map(mapfile)
    #embed(header='104 of plot')
    nside = hp.get_nside(dm_map)
    print("Nside:", nside)


    # Plot the Mollweide view of the DM map
    hp.mollview(dm_map, title=title, unit=r'pc $cm^{-3}$', min=min, max=max, cmap='Blues',)
    hp.graticule()
    plt.show()

    return None

if __name__ == '__main__':
    plot_mollwiede_view_dm_ism()
