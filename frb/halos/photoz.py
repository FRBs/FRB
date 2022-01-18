import numpy as np, os, glob

from astropy.table import Table, vstack, join
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import sigma_clipped_stats

from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator
from scipy.sparse import lil_matrix, save_npz

from frb.halos.models import ModifiedNFW, halomass_from_stellarmass
from frb.frb import FRB
from frb.galaxies import cigale as frbcig
from frb.galaxies import eazy as frb_ez
from frb.surveys import des
from frb import defs

try:
    from pathos.multiprocessing import ProcessingPool as Pool
except ImportError:
    print("You will need to run 'pip install pathos' to use some functions in this module.")
try:
    import progressbar
except ImportError:
    print("You will need to run 'pip install progressbar2' to use some functions in this module.")
try:
    from threedhst import eazyPy as ez
except ImportError:
    print("You will need to run 'pip install threedhst' to read EAZY output.")

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
                                  radius. These stars will be removed. Must contain at least two
                                  columns: "ra" and "dec". The values must be in decimal degrees
                                  and the column names are case sensitive. 
        write (bool, optional): Write output table to file?
        outfile (str, optional): Path to the output file. If not given and write is True,
                                 the table will be written to "photom_cat_J{coords}_{radius}arcmin.fits"
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
    photom_cat = cat[(cat['star_flag_r']<starflagval)&(cat['DES_r']>starbright)]

    # Remove GAIA stars if given
    if gaiacat:
        gaia_tab = Table.read(gaiacat)
        gaia_coords = SkyCoord(gaia_tab['ra'], gaia_tab['dec'], unit="deg")
        idx, d2d, _ = gaia_coords.match_to_catalog_sky(des_coords)
        matched_des = cat[idx][d2d<1*u.arcsec]
        matched_gaia = gaia_tab[d2d<1*u.arcsec]
        photom_cat = Table(np.setdiff1d(photom_cat, matched_des))

    if write:
        if outfile is None:
            coordstr = coords.to_string(style='hmsdms', sep="", precision=2).replace(" ", "")
            outfile = "photom_cat_J{:s}_{:0.1f}_arcmin.fits".format(coordstr,radius.to('arcmin').value) 
        photom_cat.write(outfile, overwrite=True)
    
    return photom_cat

def _gen_eazy_tab(photom_cat:Table, input_dir:str="eazy_in", name:str="FRB180924", out_dir:str="eazy_out", output_tab:str="no_stars_eazy.fits")->Table:
    """
    Run EAZY on the photometry and produce p(z) estimates.
    Args:
        photom_cat (Table): Photometry catalog.
        input_dir (str, optional): Folder where EAZY config files are written.
        name (str, optional): frb name. Will be passed onto frb.galaxies.eazy.eazy_input_files
            to generate input files.
        out_dir (str, optional): Folder where EAZY output is stored.
        output_tab (str, optional): Name of the output summary table fits file.
    Returns:
        joined_tab (Table): EAZY results table joined (type:inner) with photom_cat based
            on the id/ID columns. 
    """
    
    # Prepare EAZY
    frb_ez.eazy_input_files(photom_cat, input_dir, name, out_dir,
                            prior_filter="r", zmin=0.01)
    
    # Run it
    logfile = os.path.join(out_dir, "eazy_run.log")
    frb_ez.run_eazy(input_dir, name, logfile)

    # read EAZY output
    photz_file = os.path.join(out_dir, "photz.zout")
    eazy_tab = Table.read(photz_file, format="ascii")
    eazy_tab.rename_column('id','ID')
    
    # Combine the input catalog with EAZY output
    joined_tab = join(photom_cat, eazy_tab, 'ID')
    joined_tab.write(output_tab, overwrite=True)

    return joined_tab

def _create_cigale_in(photom_cat:Table, zmin:float = 0.01, zmax:float=0.35, n_z:int = 35, cigale_input:str = "cigin_minz_zfrb.fits")->Table:
    """
    Take the photometry table and
    create a new table with redshifts.
    For each galaxy, create multiple entries
    with different redshifts from 0 to 2.
    These redshifts will be uniformly spaced.
    Args:
        photom_cat (Table): Photometry catalog
        zmin (float, optional): Minimum redshift for analysis.
        zmax (float, optional): Maximum redshift for analysis.
        n_z (int, optional): Number of redshift grid points.
        cigale_input (str, optional): Name of input file to be produced.
    Returns:
        stacked_photom (Table): A table with multiple groups, one for each galaxy.
            Each entry in a group has the same photometry but different redshift values.
            This way, CIGALE can be run on the same galaxy at multiple redshift guesses
            in one go.
    """
    # Define z values
    z_range = np.linspace(zmin, zmax, n_z)

    photom_cat['redshift'] = z_range[0] # Set up initial redshift value
    photom_cat['ID'] = photom_cat['ID'].astype(str) # Convert form int to str
    photom_cat.sort("separation")
    photom_cat['ID'] = [ID.zfill(5)+"_{:0.2f}".format(z_range[0]) for ID in photom_cat['ID']]

    # Create new table
    stacked_photom = photom_cat.copy()
    for z in z_range[1:]:
        newphotom = photom_cat.copy()
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

def _gen_cigale_tab(stacked_photom:Table, n_chunks:int=10, n_cores:int=25, outdir:str=DEFAULT_DATA_FOLDER)->Table:
    """
    Run CIGALE and produce a table of results.
    Args:
        stacked_photom (Table): Table with a group for each galaxy. Output of _create_cigale_in.
        n_chunks (int, optional): How many chunks do you want to split stacked_photom.groups into?
            Just so that galaxies are not redone in case of a crash.
        n_cores (int, optional): Number of CPU threads to be used.
        outdir (str, optional): Path to the output directory.
    Returns:
        full_results (Table): CIGALE output with stellar mass and error for all entries
            in stakced_photom.
    """
    chunk_size = int(len(stacked_photom.groups)/n_chunks)

    # Only compute SFH and Stellar mass.
    compute_variables = ['stellar.m_star']

    for num in range(n_chunks):
        cigale_outdir = os.path.join(outdir,"out_minz_zfrb_chunk{}".format(num))
        # Check if a chunk has already been computed
        if os.path.isdir(cigale_outdir):
            print("Chunk {} has already been analyzed.".format(num))
            continue
        else:
            cig_photom = stacked_photom.groups[num*chunk_size:(num+1)*chunk_size]
            # Run cigale on each chunk of galaxies.
            frbcig.run(cig_photom, 'redshift', plot=False,
                outdir=cigale_outdir, cores=n_cores, variables=compute_variables, save_sed=False)
    
    # Read and combine the CIGALE results
    cigfolders = glob.glob(os.path.join(outdir, "out_minz_zfrb_chunk*"))
    relevant_cols = ['id', 'bayes.stellar.m_star', 'bayes.stellar.m_star_err']
    all_results = []
    for folder in cigfolders:
        results = Table.read(os.path.join(folder, "results.fits"))
        all_results.append(results[relevant_cols])

    full_results = vstack(all_results)
    full_results.write(os.path.join(outdir, "cigale_full_output.fits"), overwrite=True)
    return full_results 

def _load_cigale_results(cigale_input:str, cigale_output:str)->Table:
    """
    Load the CIGALE stellar mass data.
    Args:
        cigale_input (str): cigale input file path.
        cigale_output (str): cigale_output file path.
    Returns:
        trim_tab (Table): Summary table with CIGALE results.
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

def _sample_eazy_redshifts(gal_ID:int, eazy_outdir:str, ndraws:int = 1000)->np.ndarray:
    """
    Returns a sample of redshifts drawn from the
    EAZY photo-z PDF of galaxy <gal_iD>.
    Args:
        gal_ID(int): ID number of the galaxy in the EAZY table.
        eazy_outdir(str): Path to the EAZY results folder
        ndraws(int, optional): Number of redshift samples desired.
    Returns:
        sample_z (np.ndarray): Redshift sample array of length ndraws.
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

def _mhalo_lookup_table(z:float, npz_out:str = "m_halo_realizations", n_cores:int = 8):
    """
    For a given z, produce realizations of m_halo for relevant
    m_star values using only the uncertainty in the SHMR relation.
    Internal function. Use directly if you know what you're doing.
    Args:
        z (float): redshift
        npz_out(str, optional): output .npz file path.
        n_cores(int, optional): Number of CPU threads used for parallel processing.
    """

    # Define a range of stellar masses
    n_star = 1000
    log_mstar_array = np.linspace(6, 11, n_star)

    # Instantiate a 2D array
    n_halo = 10000
    log_mhalo_array = np.zeros((n_star, n_halo))

    def mhalo_factory(log_mstar:float, z:float, n_cores = n_cores)->np.ndarray:
        """
        Parallelize m_halo computations for a given log_mstar array.
        """    
        p = Pool(n_cores)
        func = lambda x: halomass_from_stellarmass(x, z = z, randomize=True)
        log_mhalo_array = p.map(func, log_mstar)
        
        return log_mhalo_array


    # Loop over log_mstar:
    for idx, log_mstar in enumerate(log_mstar_array):
        temp_log_mstar = np.full(n_halo, log_mstar)
        log_mhalo_array[idx] = mhalo_factory(temp_log_mstar, z = z, n_cores = n_cores)
    
    # Store this in an .npz file
    np.savez_compressed(npz_out, MSTAR=log_mstar_array, MHALO=log_mhalo_array)
    return

def mhalo_lookup_tables(z_grid:list, datafolder:str=DEFAULT_DATA_FOLDER, n_cores:int=8):
    """
    For each z in z_grid, produces a fits file containing m_halo values
    corresponding to a fixed grid of m_star values. The values are produced
    by sampling the Moster+13 SHMR relation. The fits files can then be
    used to produce interpolation functions of the moments of the m_halo
    distribution (e.g. mean, std.dev) as a function of redshift and log_mstar.
    Args:
        z_grid (list or np.ndarray): List of redshift values to be sampled.
        datafolder (str, optional): Path to the directory where the results will be stored.
        n_cores (int, optional): Number of CPU threads used for parallel processing.
    """

    # Just loop over z_grid and produce the fits files.
    for z in z_grid:
        realization_file = os.path.join(datafolder, "mhalo_realization_z_{:0.2f}".format(z))
        _mhalo_lookup_table(z, realization_file, n_cores)

    return

def _mhalo_realizations(log_mstar:float, log_mstar_err:float, z:float,
                        mean_interp:interp2d, stddev_interp:interp2d,
                        n_mstar:int=100, n_norm:int=10, max_log_mhalo:float=12.8)->np.ndarray:
    """
    Using the lookup tables generated (see function mhalo_lookup_tables), produce
    realiztions of mhalo. This takes into account both the stellar mass uncertainty
    and the uncertainty in the SMHR relation from Moster+13.
    Args:
        log_mstar (float): log stellar mass in M_sun.
        log_mstar_err (float): log error in log_mstar
        z (float): redshift
        mean_interp (interp2d): <log_mhalo(log_mstar, z)> (based on SHMR)
        stddev_interp (interp2d): std.dev. log_mhalo(log_mstar, z) (based on SHMR)
        n_mstrar (int, optional): Number of m_star samples to be produced.
        n_norm (int, optional): Number of m_halo samples for each m_star sample.
        max_log_mhalo (float, optional): Maximum allowed log halo mass. log halo masses
            are capped artificially to this value if any exceed.
    Returns:
        mhalo_reals (np.ndarray): log_mhalo realizations.
    """

    # First produce realizations of mstar from a normal distribution.
    mstar_reals = np.random.normal(log_mstar, log_mstar_err, n_mstar)

    # Then get mean values of halo masses for each stellar mass.
    mean_mhalo_reals = mean_interp(mstar_reals, z)
    mean_mhalo_reals = np.minimum(mean_mhalo_reals, max_log_mhalo) # Set a cutoff for the mean halo mass

    # Then get the std. dev of the halo masses for each stellar mass.
    stddev_mhalo_reals = stddev_interp(mstar_reals, z)

    # Finally, produce mhalo realizations assuming a normal distribution
    # with the means and std.devs from above.
    dummy_normal = np.random.normal(0,1, (n_norm,n_mstar))
    mhalo_reals = np.ravel(stddev_mhalo_reals*dummy_normal+mean_mhalo_reals)

    return mhalo_reals

def _dm_pdf(cigale_tab:Table, eazy_outdir:str,
            mean_interp:interp2d, stddev_interp:interp2d,
            ang_dia_interp:interp1d, dm_interpolator:RegularGridInterpolator,
            n_cores:int = 8):
    """
    For a given galaxy, compute its PDF of
    DM from the CIGALE and EAZY inputs.
    Args:
        cigale_tab (Table): On of the groups
            from the full cigale result. This
            group contains data on only one galaxy
            at various assumed redshifts. 
        eazy_outdir (str): Path to the directory with EAZY output
        mean_interp (interp2d): <log_mhalo(log_mstar, z)> (based on SHMR)
        stddev_interp (interp2d): std.dev. log_mhalo(log_mstar, z) (based on SHMR)
        ang_dia_interp (interp1d): angular_diameter_distance(z) (default Repo cosmology)
        dm_interpolator (RegularGridInterpolator): DM(z, offset_kpc, log_mhalo)
        n_cores (int, optional): Number of CPU threads to use.
    Returns:
        dm_values (np.ndarray): Array containing DM realizations for the galaxy.
        z_draws (np.ndarray): Array containing redshift draws from which dm_values were produced.
    """
    
    # Prepare interpolation functions from the
    # CIGALE table
    log_mstar_interp = interp1d(cigale_tab['redshift'], cigale_tab['log_mstar'], bounds_error=False, fill_value=1)
    log_mstar_err_interp = interp1d(cigale_tab['redshift'], cigale_tab['log_mstar_err'], bounds_error=False, fill_value=1)

    # Get 1000 random redshift draws from EAZY
    z_draws = _sample_eazy_redshifts(cigale_tab['gal_ID'][0], eazy_outdir)
    if np.isscalar(z_draws):
        return -99.

    # Convert the photo-z draws to mean stellar masses and errors
    log_mstar_array = log_mstar_interp(z_draws)
    log_mstar_err_array = log_mstar_err_interp(z_draws)

    func = lambda idx: _mhalo_realizations(log_mstar_array[idx], log_mstar_err_array[idx], z_draws[idx], mean_interp, stddev_interp)

    # Draw stellar mass values from a normal distribution and produce halo
    # masses, halo_mass errors
    p = Pool(n_cores)
    log_mhalos = p.map(func, np.arange(len(z_draws)))
    zz_draws = np.repeat(z_draws, len(log_mhalos[0]))
    offsets = ang_dia_interp(z_draws)*cigale_tab['sep_ang'][0]*u.arcmin.to('rad')
    oo_draws = np.repeat(offsets, len(log_mhalos[0]))
    dm_values = dm_interpolator((zz_draws, oo_draws, np.concatenate(log_mhalos)))

    return dm_values, z_draws.astype('float32') # Save memory by switching to a 32 bit representation.

def dm_grid(frb_z:float, n_z:int = 100, n_o:int = 100, n_m:int =100, max_log_mhalo:float=12.8,
            outdir:str=DEFAULT_DATA_FOLDER, outfile:str=None)->None:
    """
    Produce DM estimates for a 3D grid of
    redshift, offsets and log_halo_masses and write
    them to disk.
    Args:
        frb_z(float): frb redshift
        n_z(int, optional): size of the redshift grid. i.e. np.linspace(0, frb_z, n_z)
        n_o(int, optional): size of the offset grid. i.e. np.linspace(0, 600, n_o)
        n_m(int, optional):size of the log_halo_mass grid. i.e. np.linspace(8, 16, n_m)
        max_log_mhalo (float, optional): DM for halo masses larger than this are currently
            set to -99.0 to prevent weirdly large DM contributions from galactic halos. 
        outdir(str, optional): data directory to store results
        outfile(str, optional): name of results .npz file (within outdir).
    """
    # Redshift grid
    redshifts = np.linspace(0, frb_z, n_z)

    # Offset grid
    offsets = np.linspace(0, 600, n_o)

    # Mass grid
    log_halo_masses = np.linspace(8, 16, n_m)

    ZZ, OO, MM = np.meshgrid(redshifts, offsets, log_halo_masses, indexing='ij')
    raveled_z = ZZ.ravel()
    raveled_o = OO.ravel()
    raveled_m = MM.ravel()

    def halo_dm(idx):
        if raveled_m[idx] > max_log_mhalo: # Not necessary but just in case.
            return -99.0
        else:
            mnfw = ModifiedNFW(raveled_m[idx], alpha = 2, y0 = 2, z = raveled_z[idx])
        return mnfw.Ne_Rperp(raveled_o[idx]*u.kpc).to('pc/cm**3').value/(1+raveled_z[idx])

    p = Pool(8)

    raveled_dm = np.array(p.map(halo_dm, np.arange(n_z*n_o*n_m)))
    # Dm grid
    dm_grid = raveled_dm.reshape((n_z, n_o, n_m))
    if not outfile:
        outfile = os.path.join(outdir, "halo_dm_data.npz")

    np.savez_compressed(outfile, redshifts=redshifts, offsets=offsets, m_halo=log_halo_masses, dm=dm_grid)

    return

def _instantiate_intepolators(datafolder:str=DEFAULT_DATA_FOLDER, dmfilename:str=None, frb_name:str="FRB180924")->list:
    """
    Produce interpolator functions
    for key quantities required
    for the analysis.
    Args:
        datfolder(str, optional): Folder where the interpolation data files exist
        dmfilename(str, optional): file name (within datafolder) for the DM interpolation data.
        frb_name(str, optional): Assumes "FRB180924" by default.
    Returns:
        dm_interpolator (RegularGridInterpolator): DM(z, offset_kpc, log_mhalo)
        mean_interp (interp2d): <log_mhalo(log_mstar, z)> (based on SHMR)
        stddev_interp (interp2d): std.dev. log_mhalo(log_mstar, z) (based on SHMR)
        ang_dia_interp (interp1d): angular_diameter_distance(z) (default Repo cosmology)
    """

    # DM for a variety of halo parameters.
    if not dmfilename:
        dmfilename = "halo_dm_data.npz" 
    dmdata = np.load(dmfilename)
    redshifts = dmdata['redshifts']
    offsets = dmdata['offsets']
    log_mhalos = dmdata['m_halo']
    dm_grid = dmdata['dm']

    dm_interpolator = RegularGridInterpolator((redshifts, offsets, log_mhalos), dm_grid,bounds_error=False, fill_value=0.)

    # Halo mass mean and variance from stellar mass
    frb = FRB.by_name(frb_name)

    realization_files = glob.glob(os.path.join(datafolder, "mhalo_realization_z*.npz"))
    realization_files.sort()

    # Define redshift grid
    zgrid = np.linspace(0, frb.z, 10)
    
    # Now initialize arrays to store mean and std.dev.
    mean_arrays = []
    stddev_arrays = []
    # Loop through files, compute mean & std.dev of log_mhalo for log_mstar
    for file in realization_files:
        loaded = np.load(file)
        log_mhalo = loaded['MHALO']
        mean_mhalo, _, stddev_mhalo = sigma_clipped_stats(log_mhalo, sigma = 20, axis=1)
        mean_arrays.append(mean_mhalo)
        stddev_arrays.append(stddev_mhalo)
    
    # laoded is going to be from the last file in the loop. The first entry contains
    # a stellar mass array.
    log_mstar = loaded['MSTAR']
    mean_interp = interp2d(log_mstar, zgrid, np.array(mean_arrays), bounds_error=False)
    stddev_interp = interp2d(log_mstar, zgrid, np.array(stddev_arrays), bounds_error=False)

    # Angular diameter distance
    z = np.linspace(0,7, 10000)
    ang_dia_dist = defs.frb_cosmo.angular_diameter_distance(z).to('kpc').value
    ang_dia_interp = interp1d(z, ang_dia_dist, bounds_error=False, fill_value='extrapolate')

    # Return interpolators
    return dm_interpolator, mean_interp, stddev_interp, ang_dia_interp

def dm_for_all_galaxies(frb:FRB, input_catfile:str, datafolder:str,
                        n_cores:int=8, n_gals:int = None):
    """
    Produce DM estimates for all the galaxies provided by the user. Creates
    two files : "DM_halos_zdraws.npz" which contains all the redshift draws
    used for the DM realizations and "DM_halos_final.npz" which contains the
    DM realizations themselves. Each row in each of these files corresponds to
    one galaxy and each z draw corresponds to 1000 DM realizations for a galaxy.
    Args:
        frb (FRB): The FRB object of interest.
        input_catfile (str): Path to the input catalog of photometry. Assumed
            to be from DES for now.
        datafolder (str): Path to the folder in which results will be saved.
        n_cores (int, optional): Number of CPU threads to be used for computation.
        n_gals (int, optional): Limit analysis to n_gals galaxies for testing purposes.
    
    """
    # Load the input catalog
    master_cat = Table.read(input_catfile)
    # First run EAZY on that master_cat
    print("Running EAZY on the input catalog first ...")
    eazy_outdir = os.path.join(datafolder, "eazy_output")
    eazy_tab = _gen_eazy_tab(master_cat, datafolder, frb.frb_name, eazy_outdir)
    print("Done")

    # Create a CIGALE input file
    print("Creating a CIGALE input file...")
    stacked_photom = _create_cigale_in(master_cat, zmax = frb.z+0.03)
    print("Running CIGALE ...")
    cigale_output = _gen_cigale_tab(stacked_photom, outdir=datafolder, n_cores=n_cores)
    # Load CIGALE results
    cigale_input = input_catfile
    cigale_output = os.path.join(datafolder,"cigale_full_output.fits")
    cigale_tab = _load_cigale_results(cigale_input, cigale_output)
    print("CIGALE results loaded.")

    # Prepare interpolator functions
    dm_interpolator, mean_interp, stddev_interp, ang_dia_interp = _instantiate_intepolators(datafolder)
    print("Interpolators created.")

    # Reduce the sample size for testing purposes.
    if (n_gals!=None) & (type(n_gals)==int):
        eazy_tab = eazy_tab[:n_gals]
    
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
                dm_realizations[idx], z_draws[idx] = _dm_pdf(cigale_galaxy, eazy_outdir, mean_interp,
                                        stddev_interp, ang_dia_interp, dm_interpolator,
                                        n_cores = 20)
                
            bar.update(idx)
    # Save results to file
    np.savez_compressed(os.path.join(datafolder, "DM_halos_zdraws.npz"), z_draws=z_draws)
    save_npz(os.path.join(datafolder,"DM_halos_final.npz"), dm_realizations.tocsr())
    print("Done calculating")

    return

