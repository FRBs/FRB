import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from matplotlib import pyplot as plt
import pdb
from frb.surveys import des
import sys,os
import progressbar as pb #pip install progressbar2
from matplotlib import pyplot as plt
import seaborn as sns
def get_catalog(coords,size=1*u.deg):
    """
    Download a catalog objects within
    a square of input `size` centered
    around `coords`.

    Args:
        coords (astropy SkyCoord): central coordinates
        size (astropy Angle, optional): Size of the square FoV around
                              the central coordinates
    Returns:
        catalog (astropy Table): DES DR1 search results
    """
    survey = des.DES_Survey(coords,size/np.sqrt(2))
    catalog =  survey.get_catalog(print_query=False)
    select = (catalog['ra']>coords.ra.value-size.value/2)&(catalog['ra']<coords.ra.value+size.value/2)
    select = select*(catalog['dec']>coords.dec.value-size.value/2)&(catalog['dec']<coords.dec.value+size.value/2)
    catalog = catalog[select]
    return catalog

def _generate_coord_grid(coords,size=1*u.deg,resolution=3600):
    """
    Genereate a unifrom grid
    of SkyCoords centered around `coords` within
    an area `size`x`size` with a default `resolution`
    of 1000 points along each axis.

    Args:
        coords (astropy SkyCoord): central coordinates
        size (astropy Angle, optional): Size of the square FoV around

    Returns:
            SkyCoord:
    """
    ra,dec = coords.ra.value,coords.dec.value
    ra_arr = np.linspace(ra-size.value/2,ra+size.value/2,resolution)
    dec_arr = np.linspace(dec-size.value/2,dec+size.value/2,resolution)
    rr,dd = np.meshgrid(ra_arr,dec_arr)
    return SkyCoord(rr.ravel(),dd.ravel(),unit="deg")

def get_frac_within_sep(coords,catalog,sep=1*u.arcsec,resolution=1000,size=1*u.deg,band='r',crit_mag=22):
    """
    Obtain the fraction of sightlines on the unifrom grid
    defined in _generate_coord_grid falling within `sep` distance
    of a galaxy in `catalog` with mag less than `crit_mag` in `band`.
    The catalog and grid are defined in a square centered around `coords`
    of side length `size` and grid has input linear `resolution`. 
    """
    grid = _generate_coord_grid(coords,size,resolution)
    catalogcoord = SkyCoord(catalog['ra'],catalog['dec'],unit="deg")
    idx, sep2d, _ = grid.match_to_catalog_sky(catalogcoord)
    newtab = Table()
    newtab['gridcoord'] = grid
    newtab['sep'] = sep2d
    newtab["DES_"+band] = catalog[idx]["DES_"+band]
    #compute frac
    select = (newtab["DES_"+band]<crit_mag)&(newtab['sep']<sep.to(u.deg).value)
    frac = np.sum(select)/len(grid)
    return frac

def random_sightlines(n=100,resolution=3600,sep=1*u.arcsec,size=1*u.deg,
                        band='r',crit_mag=22,outfile="random_sights.txt"):
    """
    Query a contigous quare patch of DES `n` times to obtain output
    from `get_frac_within_sep` and store it to `outfile`.
    """
    ra = 22.5 + np.random.rand(n)*45 #limit RA to [22.5deg,67.5deg]
    dec = -60 + np.random.rand(n)*30 #limit DEC to [-60deg,-30deg]
    rand_coords = SkyCoord(ra,dec,unit="deg")
    fracs = np.zeros(n)
    bar = pb.ProgressBar(max_value=n)
    bar.start()
    for num,coords in enumerate(rand_coords):
        catalog = get_catalog(coords,size)
        fracs[num] = get_frac_within_sep(coords,catalog,sep=sep,resolution=resolution,band=band,crit_mag=crit_mag)
        bar.update(num+1)
    #Save to file
    np.savetxt(outfile,fracs)
    return fracs
def plot_hist(fracs,bins=5):
    """
    Plot a histogram of fractions obtained from `random_sightlines`
    """
    sns.set(font_scale=1.3,font="serif",style="ticks")
    sns.axes_style(rc={"weight":"bold"})
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,8))
    ax.set_xlabel("Percentage of sightlines within 1'' of a galaxy (r<22)",labelpad=20,weight="bold")
    ax.set_ylabel("Number of DES patches (1 sq. deg)",labelpad=20,weight="bold")
    ax.hist(fracs*100,bins=bins,edgecolor="k",color="#366293")
    median_frac = np.median(fracs)*100 
    ax.axvline(x=median_frac,linestyle="--",color="#262626")
    ax.annotate("Median = {:.3f}".format(median_frac),(median_frac+0.01,22),color="#262626")
    plt.savefig("/home/sunil/Desktop/DES_fracs.png",dpi=300,bbox_inches="tight",pad_inches=0.1)
    plt.show()


    
"""
If you're running this for the first time or you 
want to regenerate the database, uncomment the following line
and run.
"""
#fracs = random_sightlines(n=100)

fracs = np.loadtxt("random_sights.txt")
plot_hist(fracs,bins=9)


