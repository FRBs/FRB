import numpy as np, os, glob

from astropy.io import fits
from astropy.table import Table, vstack, join
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import Planck15 as p15
from astropy import visualization as vis
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats

from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator
from scipy.sparse import lil_matrix, save_npz
from scipy import stats as sci_st
from scipy.integrate import simps

from frb.halos.models import ModifiedNFW, halomass_from_stellarmass
from frb.frb import FRB
from frb.galaxies import cigale as frbcig
from frb.galaxies import eazy as frb_ez
from frb.surveys import des
from frb.surveys.catalog_utils import convert_mags_to_flux

from pathos.multiprocessing import ProcessingPool as Pool
import progressbar

from matplotlib import pyplot as plt

from threedhst import eazyPy as ez

from importlib import reload

DEFAULT_DATA_FOLDER = "data"

def get_des_data(coords:SkyCoord, radius:u.Quantity=15.*u.arcmin, starbright:float=17,
                 starflagval:float=0.9, gaiacat:str=None, write:bool=False, outfile:str=None)->Table:
    """
    Download photometry for galaxies within an FRB field.
    Args:
        coords (SkyCoord): Coordinates of the center of a cone search.
        radius (Quantity, optional): Radius of cone search.
        starbright (float, optional): Lower limit of r band mag. Objects brighter
                                    than this will be removed.
        starflagval (float, optional): Upper limit for a morphology-based classifier
                                       flag. Objects more point-like (i.e. higher value)
                                       will be filtered out.
        gaicat (str, optional): Optional file with gaia catalog of stars within the same search
                                  radius. These stars will be removed.
        write (bool, optional): Write output table to file?
        outfile (str, optional): Path to the output file. If not given and write is True,
                                 the table will be written to "DES_cat_J{coords}_{radius}arcmin.fits"
                                 in the current working directory.
    Returns:
        des_data (Table): Table of DES galaxies within the search radius.
    """
    # Download catalog
    survey = des.DES_Survey(coords, radius)
    cat = survey.get_catalog()

    # Add separation info
    des_coords = SkyCoord(cat['ra'],cat['dec'], unit="deg")
    dessep = coords.separation(des_coords).to('arcmin').value
    cat['separation'] = dessep
    cat.sort("separation")
    cat_colnames = cat.colnames

    # Add a convenient unique ID
    cat['ID'] = np.arange(len(cat))+1
    cat = cat[['ID']+cat_colnames]

    # Make brightness and morphology cuts
    des_cat = cat[(cat['star_flag_r']<starflagval)&(cat['DES_r']>starbright)]

    # Remove GAIA stars if given
    if gaiacat:
        gaia_tab = Table.read(gaiacat)
        gaia_coords = SkyCoord(gaia_tab['ra'], gaia_tab['dec'], unit="deg")
        idx, d2d, _ = gaia_coords.match_to_catalog_sky(des_coords)
        matched_des = cat[idx][d2d<1*u.arcsec]
        matched_gaia = gaia_tab[d2d<1*u.arcsec]
        des_cat = Table(np.setdiff1d(des_cat, matched_des))

    if write:
        if outfile is None:
            coordstr = coords.to_string(style='hmsdms', sep="", precision=2).replace(" ", "")
            outfile = "DES_cat_J{:s}_{:0.1f}_arcmin.fits".format(coordstr,radius.to('arcmin').value) 
        des_cat.write(outfile, overwrite=True)
    
    return des_cat

def gen_eazy_tab(des_cat:Table, input_dir:str="eazy_in", name:str="FRB180924", out_dir="eazy_out", output_tab="no_stars_eazy.fits"):
    """
    Run EAZY on the photometry and produce p(z) estimates.
    """
    
    # Prepare EAZY
    frb_ez.eazy_input_files(des_cat, input_dir, name, out_dir,
                            prior_filter="r", zmin=0.01)
    
    # Run it
    logfile = os.path.join(out_dir, "eazy_run.log")
    frb_ez.run_eazy(input_dir, name, logfile)

    # read EAZY output
    photz_file = os.path.join(out_dir, "photz.zout")
    eazy_tab = Table.read(photz_file, format="ascii")
    eazy_tab.rename_column('id','ID')
    
    # Combine the input catalog with EAZY output
    joined_tab = join(des_cat, eazy_tab, 'ID')
    joined_tab.write(output_tab, overwrite=True)

    return joined_tab

def create_cigale_in(des_cat:Table, zmin:float = 0.01, zmax:float=0.35, n_z:int = 35, cigale_input:str = "cigin_minz_zfrb.fits"):
    """
    Take the photometry table and
    create a new table with redshifts.
    For each galaxy, create multiple entries
    with different redshifts from 0 to 2.
    These redshifts will be uniformly spaced.
    """
    # Define z values
    z_range = np.linspace(zmin, zmax, n_z) # Go up in steps of 0.01

    des_cat['redshift'] = z_range[0] # Set up initial redshift value
    des_cat['ID'] = des_cat['ID'].astype(str) # Convert form int to str
    des_cat.sort("separation")
    des_cat['ID'] = [ID.zfill(5)+"_{:0.2f}".format(z_range[0]) for ID in des_cat['ID']]

    # Create new table
    stacked_photom = des_cat.copy()
    for z in z_range[1:]:
        newphotom = des_cat.copy()
        newphotom['redshift'] = z
        for entry in newphotom:
            entry['ID'] = entry['ID'].replace("_0.01", "_{:0.2f}".format(z))
        stacked_photom = vstack([stacked_photom, newphotom])
    
    
    # Sort table by ID
    stacked_photom = stacked_photom.group_by('ID')

    # Write to disk
    stacked_photom.write(cigale_input, overwrite=True)
    print("Wrote to disk {:s}".format(cigale_input))
    return stacked_photom

def gen_cigale_tab(stacked_photom:Table, n_groups:int=10, n_cores:int=25):
    chunk_size = int(len(stacked_photom.groups)/n_groups)

    # Only compute SFH and Stellar mass.
    compute_variables = ['sfh.sfr','sfh.tau_main','sfh.age_main','stellar.m_star']

    for num in range(n_groups):
        cigale_outdir = os.path.join(DEFAULT_DATA_FOLDER,"out_minz_zfrb_chunk{}".format(num))
        # Check if a chunk has already been computed
        if os.path.isdir(cigale_outdir):
            print("Chunk {} has already been analyzed.".format(num))
            continue
        else:
            cig_photom = stacked_photom.groups[num*chunk_size:(num+1)*chunk_size]
            # Run cigale on each chunk of galaxies.
            frbcig.run(cig_photom, 'redshift', plot=False,
                outdir=cigale_outdir, cores=n_cores, variables=compute_variables, save_sed=False)
    return

def create_interpolators(datafolder:str, frb_name:str="FRB180924"):
    """
    Produce interpolator functions
    for key quantities required
    for the analysis.
    """

    # DM for a variety of halo parameters.
    hdulist = fits.open(os.path.join(datafolder, "halo_dm_data.fits"))
    redshifts = hdulist[1].data
    offsets = hdulist[2].data
    log_mhalos = hdulist[3].data
    dm_grid = hdulist[4].data

    dm_interpolator = RegularGridInterpolator((redshifts, offsets, log_mhalos), dm_grid,bounds_error=False, fill_value=0.)

    # Halo mass mean and variance from stellar mass
    frb = FRB.by_name(frb_name)

    realization_files = glob.glob(os.path.join(datafolder, "mhalo_realization_z*.fits"))
    realization_files.sort()

    # Define redshift grid
    zgrid = np.linspace(0, frb.z, 10)
    
    # Now initialize arrays to store mean and std.dev.
    mean_arrays = []
    stddev_arrays = []
    # Loop through files, compute mean & std.dev of log_mhalo for log_mstar
    for file in realization_files:
        hdulist = fits.open(file)
        log_mhalo = hdulist[2].data
        mean_mhalo, _, stddev_mhalo = sigma_clipped_stats(log_mhalo, sigma = 20, axis=1)
        mean_arrays.append(mean_mhalo)
        stddev_arrays.append(stddev_mhalo)
    
    # hdulist is going to be from the last file in the loop. The first HDU contains
    # stellar mass array.
    log_mstar = hdulist[1].data
    mean_interp = interp2d(log_mstar, zgrid, np.array(mean_arrays), bounds_error=False)
    stddev_interp = interp2d(log_mstar, zgrid, np.array(stddev_arrays), bounds_error=False)

    # Angular diameter distance
    z = np.linspace(0,7, 10000)
    ang_dia_dist = p15.angular_diameter_distance(z).to('kpc').value
    ang_dia_interp = interp1d(z, ang_dia_dist, bounds_error=False, fill_value='extrapolate')

    # Return interpolators
    return dm_interpolator, mean_interp, stddev_interp, ang_dia_interp

def load_cigale_results(cigale_input, cigale_output):
    """
    Load the CIGALE stellar mass data.
    """
    cigin = Table.read(cigale_input)
    cigtab = Table.read(cigale_output)
    
    # Trim the output table
    trim_tab = cigtab[['id', 'bayes.stellar.m_star', 'bayes.stellar.m_star_err']]

    # produce some extra columns
    trim_tab['redshift'] = 0.0
    trim_tab['gal_ID'] = 1

    for entry in trim_tab:
        entry['gal_ID'] = int(entry['id'][:-5])
        entry['redshift'] = float(entry['id'][-4:])

    # produce a column for angular separation
    trim_tab.sort('id')
    trim_tab = trim_tab.group_by('gal_ID')
    trim_tab['sep_ang'] = 99.0

    for group in trim_tab.groups:
        group['sep_ang'] = cigin['separation'][cigin['ID'] == group['gal_ID'][0]][0]

    # A similar column for separation in kpc
    #trim_tab['sep_kpc'] = p15.angular_diameter_distance(trim_tab['redshift']).to('kpc').value*trim_tab['sep_ang']*u.arcmin.to('rad')

    # Rename the stellar mass columns
    trim_tab.rename_columns(['bayes.stellar.m_star', 'bayes.stellar.m_star_err'],['log_mstar', 'log_mstar_err'])
    # Convert to logarithmic values
    trim_tab['log_mstar_err'] = (np.log10(trim_tab['log_mstar']+trim_tab['log_mstar_err']) - 
                                np.log10(np.abs(trim_tab['log_mstar']-trim_tab['log_mstar_err'])))/2
    trim_tab['log_mstar'] = np.log10(trim_tab['log_mstar'])

    return trim_tab

def sample_eazy_redshifts(gal_ID, eazy_outdir, ndraws = 1000):
    """
    Returns a sample of redshifts drawn from the
    EAZY photo-z PDF of galaxy <gal_iD>.
    """
    # Get posterior
    zgrid, pz = ez.getEazyPz(gal_ID-1,OUTPUT_DIRECTORY=eazy_outdir)
    # Force a value of 0 at z = 0
    zgrid = np.hstack([[0],zgrid])
    pz = np.hstack([[0],pz])
    if np.all(np.diff(zgrid) == 0):
        return -99

    # make a CDF
    cdf_z = np.cumsum(pz)
    cdf_z /= np.max(cdf_z)
    cdf_interp = interp1d(cdf_z, zgrid, kind="linear", fill_value=0, bounds_error=False)

    # Use uniform distribution to produce random draws from the CDF
    sample_u = np.random.rand(ndraws)
    sample_z = cdf_interp(sample_u)
    return sample_z

def mhalo_realizations(log_mstar, log_mstar_err, z, mean_interp, stddev_interp, n_mstar=100, n_norm=10):
    mstar_reals = np.random.normal(log_mstar, log_mstar_err, n_mstar)
    mean_mhalo_reals = mean_interp(mstar_reals, z)
    # Set a cutoff for the mean halo mass
    mean_mhalo_reals[mean_mhalo_reals>12.8] = 12.8
    stddev_mhalo_reals = stddev_interp(mstar_reals, z)

    dummy_normal = np.random.normal(0,1, (n_norm,n_mstar))
    mhalo_reals = np.ravel(stddev_mhalo_reals*dummy_normal+mean_mhalo_reals)
    return mhalo_reals

def dm_pdf(cigale_tab, eazy_outdir, mean_interp, stddev_interp, ang_dia_interp, dm_interpolator, ncores = 8):
    """
    For a given galaxy, compute its PDF of
    DM from the CIGALE and EAZY inputs.
    Args:
        cigale_tab (Table): On of the groups
            from the full cigale result. This
            group contains data on only one galaxy
            at various assumed redshifts. 
    """
    
    # Prepare interpolation functions from the
    # CIGALE table
    log_mstar_interp = interp1d(cigale_tab['redshift'], cigale_tab['log_mstar'], bounds_error=False, fill_value=1)
    log_mstar_err_interp = interp1d(cigale_tab['redshift'], cigale_tab['log_mstar_err'], bounds_error=False, fill_value=1)

    # Get 1000 random redshift draws from EAZY
    z_draws = sample_eazy_redshifts(cigale_tab['gal_ID'][0], eazy_outdir)
    if np.isscalar(z_draws):
        return -99.

    # Convert the photo-z draws to mean stellar masses and errors
    log_mstar_array = log_mstar_interp(z_draws)
    log_mstar_err_array = log_mstar_err_interp(z_draws)

    func = lambda idx: mhalo_realizations(log_mstar_array[idx], log_mstar_err_array[idx], z_draws[idx], mean_interp, stddev_interp)



    # Draw stellar mass values from a normal distribution and produce halo
    # masses, halo_mass errors
    p = Pool(ncores)
    log_mhalos = p.map(func, np.arange(len(z_draws)))
    zz_draws = np.repeat(z_draws, len(log_mhalos[0]))
    offsets = ang_dia_interp(z_draws)*cigale_tab['sep_ang'][0]*u.arcmin.to('rad')
    oo_draws = np.repeat(offsets, len(log_mhalos[0]))
    dm_values = dm_interpolator((zz_draws, oo_draws, np.concatenate(log_mhalos)))

    return dm_values, z_draws.astype('float32')

def main(frb:FRB, datafolder:str, master_cat:Table=None, ngals:int = None):
    """
    Run the analysis on all galaxies.
    """
    

    if master_cat is None:
        # Create data file
        master_cat = get_des_data(frb.coord)
        print("Downloaded the catalog.")

    # Create a CIGALE input file
    stacked_photom = create_cigale_in(master_cat)
    print("Created a CIGALE input file.")

    # Prepare interpolator functions
    dm_interpolator, mean_interp, stddev_interp, ang_dia_interp = create_interpolators(datafolder)
    print("Interpolators created.")

    # Load CIGALE results
    cigale_input = os.path.join(datafolder, "180924_fg_DES_WISE.fits")
    cigale_output = os.path.join(datafolder,"cigresults_180924fg_DES_WISE_bayes.fits")
    cigale_tab = load_cigale_results(cigale_input, cigale_output)
    print("CIGALE results loaded.")

    # Load EAZY results
    eazy_outdir = os.path.join(datafolder, "eazy_output")
    eazy_tab = Table.read(os.path.join(datafolder, "no_stars_eazy.fits"))

    # Reduce the sample size for testing purposes.
    if (ngals!=None) & (type(ngals)==int):
        eazy_tab = eazy_tab[:ngals]
    
    # Loop through galaxies
    print("Computing DM realizations for all galaxies ...")
    # Initialize storage for the DM realizations and the redshifts at which these are computed.
    dm_realizations = lil_matrix((len(eazy_tab), 1000000))
    z_draws = np.zeros((len(eazy_tab),1000), dtype='float32')

    # Begin calculating
    with progressbar.ProgressBar(max_value=len(eazy_tab)-1) as bar:
        for idx, ez_entry in enumerate(eazy_tab):
            cigale_galaxy = cigale_tab[cigale_tab['gal_ID']==ez_entry['ID']]
            if np.any(np.isnan(cigale_galaxy['log_mstar'])):
                continue
            else:
                dm_realizations[idx], z_draws[idx] = dm_pdf(cigale_galaxy, eazy_outdir, mean_interp,
                                        stddev_interp, ang_dia_interp, dm_interpolator,
                                        ncores = 20)
                
            bar.update(idx)
    # Save results to file
    np.savez_compressed(os.path.join(datafolder, "DM_halos_zdraws.npz"), z_draws=z_draws)
    save_npz(os.path.join(datafolder,"DM_halos_final.npz"), dm_realizations.tocsr())
    #prihdu = fits.PrimaryHDU()
    #dmhdu = fits.ImageHDU(dm_realizations, name="DM_halo")
    #gal_idhdu = fits.ImageHDU(eazy_tab['ID'].data, name="Gal_ID")
    #hdulist = fits.HDUList([prihdu, gal_idhdu, dmhdu])
    #hdulist.writeto(os.path.join(datafolder,"DM_halos_final.fits"), overwrite=True)
    print("Done calculating")
