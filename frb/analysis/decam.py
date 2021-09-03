"""
Module for extracting sources and performing photometry
with DECam images
"""

from astropy import units
import numpy as np, os, glob

from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.wcs import WCS, utils as wcsutils
from astropy.table import Table, join, join_skycoord, vstack, hstack, setdiff
from astropy.coordinates import SkyCoord
from astropy import units as u

from astroquery.sdss import SDSS

import sep
from photutils import Background2D
from photutils.background import MADStdBackgroundRMS, SExtractorBackground
from photutils.segmentation import SourceCatalog, SegmentationImage

_DEFAULT_PHOTOBJ_FIELDS = ['ra', 'dec', 'objid','type']
_DEFAULT_PHOTOBJ_FIELDS += ['psfMag_'+band for band in ['u', 'g', 'r', 'i', 'z']]

def make_masks(dqmimg:np.ndarray, expimg:np.ndarray, combine:bool=False)->list:
    """
    Create masks of pixels to be ignored for the analysis.
    Args:
        dqmimg (np.ndarray): Data quality map for the image.
        expimg (np.ndarray): Exposure map for the image.
        combine (bool, optional): Merge the two masks using OR?
    Returns:
        dqm_mask (np.ndarray): Mask with pixels of DQ>0 set to True.
           These are either not observed or saturated.
        exp_mask (np.ndarray): Mask with exptime = 0 set to True.
            Regions of no exposure are masked. 
    """

    # Read
    dqm_mask = dqmimg>0
    exp_mask = expimg==0
    if not combine:
        return dqm_mask, exp_mask
    else:
        return dqm_mask+exp_mask

def prepare_bkg(img:np.ndarray, sigma_clip_level:float=3.,
                box_size:tuple=(32,32), filter_size:tuple=(3,3),
                mask:np.ndarray=None, coverage_mask:np.ndarray=None)->list:
    """
    Create local background and noise maps using the Background2D
    estimators from photutils
    Args:
        img (np.ndaray): 2D image array.
        sigma_clip_level (float, optional): Threshold in number of std.devs
            above the noise floor to clip data. Ensures bright objects
            do not contribute to the background estimate.
        box_size (tuple of ints, optional): Size of the box that will be tiled
            over the image to estimate local background statistics. See
            `photutils.background.Background2D` documentation. 
        filter_size (tuple of ints, optional): The window size of the 2D median
            filter to apply to the low-resolution background map. See
            `photutils.background.Background2D` documentation.
        mask (np.ndarray, optional): 2D boolean array of the same size as the image
            with masked pixels marked True. These pixels will not be used to
            estimate the background but it will be interpolated over them.
        coverage_mask (np.ndarray, optional): 2D boolean array with True values implying no
            exposure. Background estimates for this region will be set to 0.
    Returns:
        bkg_img (np.ndarray): Local background map using SExtractorBackground (2.5*median-1.5*mean).
        noise_floor_img (np.ndarray): Local noise floor map using MADStdBackgroundRMS.
    """

    #Prepare
    sigma_clip = SigmaClip(sigma=sigma_clip_level)
    bkg_estimator = SExtractorBackground(sigma_clip=sigma_clip)
    noise_floor_estimator = MADStdBackgroundRMS(sigma_clip=sigma_clip)

    # Estimate
    bkg = Background2D(img, box_size, filter_size=filter_size,
                    sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                    mask=mask, coverage_mask=coverage_mask)
    
    noise = Background2D(img, box_size, filter_size=filter_size,
                    sigma_clip=sigma_clip, bkg_estimator=noise_floor_estimator,
                    mask=mask, coverage_mask=coverage_mask)
    
    return bkg.background, noise.background

def _source_table(data, segm, bkg_img=None, error_img=None, **sourcecat_kwargs):
    """
    Helper function for exxtract sources
    """
    # Do photometry and produce a source catalog
    segmap = SegmentationImage(segm)
    source_cat = SourceCatalog(data, segmap, error=error_img,
                                background=bkg_img, **sourcecat_kwargs).to_table()

    return source_cat, segmap


def extract_sources(img:np.ndarray, dqmimg:np.ndarray, expimg:np.ndarray, wtmap:np.ndarray,
                    threshold:float=3.,minarea:int=5,deblend_cont:float=0.001, clean_param:float=2.,
                    bkg_kwargs:dict={}, sourcecat_kwargs:dict={})->list:
    
    """
    Use sep.extract to produce a segmentation image and subsequently produce a
    source catalog with apareture photometry performed.

    Args:
        img (np.ndarray): Image without any background subtraction.
        dqmimg (np.ndarray): Data quality map.
        expimg (np.ndarray): Exposure map.
        wtmap (np.ndarray): Weight map (i.e. inverse variance map)
        threshold (float, optional): Number of std dev above the noise floor
            which counts as a valid detection.
        minarea (int, optional): Minimum number of conjoined pixels that
            constitute a detection.
        deblend_cont (float, optional): contrast level used for object deblending.
            See sep.extract documentation.
        clean_param (float, optional): sep.extract clean_param for cleaning detections.
        bkg_kwargs (dict, optional): Additional keyword args that will be passed onto
            prepare_bkg.
        sourcecat_kwargs (dict:, optional): Additional keyword args that will be passed
            onto `photutils.segementation.SourceCatalog` for producing the catalog.
    Returns:
        source_cat (Table): Catalog of detected sources with photometry.
        segmap (SegmentaionImage): Image segmentation map
    """
    # Make image masks
    dqm_mask, _ = make_masks(dqmimg, expimg)

    # Prepare bkg images
    bkg_img, noise_img = prepare_bkg(img, mask=dqm_mask, coverage_mask=None, **bkg_kwargs)

    # Subtract background
    data_sub = img - bkg_img

    # Produce deblended segmentaiton map
    # I've decided against using he coverage masks and dq masks to extract. This produces
    # annoying concave segments which are hard to filter out.
    sep.set_extract_pixstack(5000000)
    _, segm = sep.extract(data_sub, threshold, minarea=minarea,deblend_cont=deblend_cont,
                          clean_param=clean_param, err=np.ascontiguousarray(noise_img),
                          mask=None, segmentation_map=True)
    # Do photometry and produce a source catalog
    source_cat, segmap = _source_table(data_sub,segm,bkg_img,1/np.sqrt(wtmap),**sourcecat_kwargs)

    # Clean source cat
    select = ~np.isnan(source_cat['xcentroid']) # No xcentroid
    select *= ~np.isnan(source_cat['ycentroid']) # No ycentroid
    select *= ~np.isnan(source_cat['kron_flux']) # No kron flux

    trimmed_cat = source_cat[select]
    badlabels = np.setdiff1d(source_cat['label'], trimmed_cat['label'])

    segmap.remove_labels(badlabels)
    return trimmed_cat, segmap

def process_image_file(img_file:str, dqm_file:str, exp_file:str, wt_file:str,
                       outdir:str="output", extract_sources_kwargs:dict={})->None:
    """
    Run the extraction steps on the image files and produce 
    an output folder with a source catalog and segmentation map.

    Args:
        img_file (str): Path to the imagefile
        dqm_file (str): Path to the data quality map
        exp_file (str): Path to the expsure map.
        wt_file (str): Path to the inverse variance maps
        outdir (str, optional): Path to the output directory.
        extract_sources_kwargs (dict): Additional arguments for extract_sources.
    """

    # Read files
    imghdus = fits.open(img_file, memmap=True)
    dqmhdus = fits.open(dqm_file, memmap=True)
    exphdus = fits.open(exp_file, memmap=True)
    wthdus = fits.open(wt_file, memmap=True)

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Loop over all tiles
    mega_source_cat = Table()
    for idx in np.arange(1,len(imghdus)):
        tilename = imghdus[idx].header['EXTNAME']
        print("Processing {}".format(tilename))
        img = imghdus[idx].data
        dqmimg = dqmhdus[idx].data
        expimg = exphdus[idx].data
        wtimg = wthdus[idx].data

        wcs = WCS(imghdus[idx].header)
        if 'sourcecat_kwargs' not in extract_sources_kwargs.keys():
            extract_sources_kwargs['sourcecat_kwargs']= {"wcs":wcs}
        elif 'wcs' not in extract_sources_kwargs['sourcecat_kwargs'].keys():
            extract_sources_kwargs['sourcecat_kwargs']['wcs'] = wcs

        # Extract
        source_cat, segmap = extract_sources(img,dqmimg=dqmimg,expimg=expimg,wtmap=wtimg,**extract_sources_kwargs)
        source_cat['ra'] = source_cat['sky_centroid'].ra.value
        source_cat['dec'] = source_cat['sky_centroid'].dec.value
        source_cat['tile'] = tilename
        source_cat.remove_column('sky_centroid') #Just inconvenient to use in the long run. RA, Dec better separately.

        
        if len(mega_source_cat) == 0:
            mega_source_cat = source_cat
        else:
            mega_source_cat = vstack([mega_source_cat, source_cat])

        np.savez_compressed(os.path.join(outdir, tilename+"_segmap.npz"), segmap=segmap)
    # Write
    mega_source_cat.write(os.path.join(outdir, img_file.split("/")[-1].split(".")[0]+"_source_cat.fits"), overwrite=True) # Fits for better compressibility
    print("Done!")
    return

def batch_run(input_folder:str, output_folder:str=None, extract_sources_kwargs:dict={})->None:
    """
    Run process_image_file on a batch of files inside an input folder.
    Args:
        input_folder (str): Path to the input folder with image files in the top
            level. Also must have dqm mask file, exp map and weight maps.
        output_folder (str, optional): Path to the output folder. Will dump
            results into the inputfolde rif not specified.
        extract_sources_kwargs (dict, optional): Additional arguments to
            pass onto extract_sources. 
    """

    # Check if input folder exists
    assert os.path.isdir(input_folder), "Invalid input folder path {:s}".format(input_folder)

    # Create the output folder if it doesn't exist
    if output_folder is None:
        output_folder = input_folder
    elif not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    # Check if it has the right files and populate lists with file names
    imgfiles = glob.glob(os.path.join(input_folder, "c4d_*_osj_*"))
    imgfiles.sort()

    auxfile_dict = {"wtfile":["_osw_", []], "expfile":["_ose_", []], "dqmfile":["_osd_", []]}

    for file in imgfiles:
        # Check if auxiliary files are present
        for key in auxfile_dict.keys():
            filename = file.replace("_osj_", auxfile_dict[key][0])
            assert os.path.isfile(filename), "Could not find {:s} file {:s}. Did you move or rename it?".format(key, filename)
            auxfile_dict[key][1].append(filename)
        
    
    # Process it
    for idx in range(len(imgfiles)):
        temp = imgfiles[idx].split("/")[-1].split("_osj_")
        extract_dir = os.path.join(output_folder,temp[0]+"_"+temp[1][0])
        if os.path.isdir(extract_dir):
            # Don't repeat extraction if it has been done already.
            print("Skipping "+imgfiles[idx])
            continue
        print("Processing "+imgfiles[idx])
        process_image_file(imgfiles[idx], auxfile_dict['dqmfile'][1][idx],
                             auxfile_dict['expfile'][1][idx], auxfile_dict['wtfile'][1][idx],
                             extract_dir, extract_sources_kwargs)

    print("Completed processing all files in {}".format(input_folder))
    return

def _SDSS_query(coord:SkyCoord, radius:units.Quantity,
                timeout:float=120, photo_obj_fields:list=_DEFAULT_PHOTOBJ_FIELDS)->Table:
    """
    Run a simple SDSS PhotoObj query and return a table of 
    model mags for stars.
    """
    sdss_cat = SDSS.query_region(coord, radius=radius,
                                timeout=timeout, photoobj_fields=photo_obj_fields)
    
    sdss_cat = sdss_cat[sdss_cat['type']==3]
    sdss_cat['sky_centroid'] = SkyCoord(sdss_cat['ra'],sdss_cat['dec'], unit="deg")
    return sdss_cat

def photo_zpt_run(calib_file:str, band:str="r", brightest:int=50,
                  save_tab:str=None, verbose:bool=True)->float:
    """
    Run a crude segmentation map based photometry on
    the brightest stars in a standard star exposure
    and return delta = ZPT + AIRMASS_TERM*AIRMASS
    after comparing instrument fluxes against SDSS.
    """
    # Read in file
    stdhdus = fits.open(calib_file)
    # Instantiate results table
    super_merged = Table()

    # Loop over all CCDs
    try:
        import progressbar
        bar =  progressbar.ProgressBar(max_value=len(stdhdus)-2)
        pbexists = True
    except ImportError:
        pbexists=False
    print("Processing "+calib_file+"...")
    for num,hdu in enumerate(stdhdus[1:]):

        # Prepare
        data, hdr = hdu.data, hdu.header
        wcs = WCS(hdu.header)
        center = wcsutils.pixel_to_skycoord(data.shape[1]/2, data.shape[0]/2, wcs)

        bkg_img, noise_img = prepare_bkg(data)
        data_sub = data - bkg_img

        # Extract the brightest stars
        _, segmap = sep.extract(data_sub, 3., minarea=5, deblend_cont=0.005,
                clean_param=2., err=np.ascontiguousarray(noise_img),segmentation_map=True)
        sources, _ = _source_table(data_sub, segmap, bkg_img, noise_img, wcs = wcs)
        sources.sort("segment_flux")
        sources = sources[-brightest:]

        # Run an SDSS query in the same region
        try:
            sdss_cat = _SDSS_query(center, 10*u.arcmin)
        except:
            if pbexists&verbose:
                bar.update(num)
            elif verbose:
                print("Skipped "+hdr['extname']+" because of a query failure")
            continue
        # Select stars only
        sdss_cat = sdss_cat[sdss_cat['type']==3]
        sdss_cat['sky_centroid'] = SkyCoord(sdss_cat['ra'],sdss_cat['dec'], unit="deg")

        # Cross match sources against SDSS
        join_funcs={'sky_centroid': join_skycoord(0.1 * u.arcsec)}
        merged = join(sources, sdss_cat, keys="sky_centroid", join_funcs=join_funcs)
        merged.remove_columns(['sky_centroid_1', 'sky_centroid_2'])

        # Combine results of cross-matching into one giant table 
        if len(super_merged)==0:
            super_merged = merged
        else:
            super_merged = vstack([super_merged, merged])
        
        # Need a status update?
        if pbexists&verbose:
            bar.update(num)
        elif verbose:
            print("Done processing "+hdr['extname'])

    # Write table to disk?
    if save_tab is None:
        save_tab = calib_file.split("/")[-1].split(".")[0]+"_std_results.txt"
        super_merged.write(save_tab, format="ascii.fixed_width", overwrite=True)
        if verbose:
            print("Wrote results to "+save_tab)
    # Compute delta
    delta = np.median(-2.5*np.log10(super_merged['segment_flux']/stdhdus[0].header['exptime'])-
                       super_merged['psfMag_'+band.lower()])
    return delta

def get_zpt(calib_file_1:str, calib_file_2:str, band:str, verbose:bool=True):
    """
    Run `photo_zpt_run` on two standard files
    """
    delta_1 = photo_zpt_run(calib_file_1, band=band, verbose=verbose)
    delta_2 = photo_zpt_run(calib_file_2, band=band, verbose=verbose)
    
    airmass1 = fits.getheader(calib_file_1, hdu=0)['airmass']
    airmass2 = fits.getheader(calib_file_2, hdu=0)['airmass']
    
    airmass_term = (delta_1-delta_2)/(airmass1-airmass2)
    
    zpt = delta_1-airmass_term*airmass1
    
    return zpt, airmass_term

def _custom_match_func(coord1:SkyCoord, coord2:SkyCoord, seplimit:units.Quantity)->tuple:
    idx, d2d, d3d = coord1.match_to_catalog_sky(coord2)
    match = d2d<seplimit
    idx1 = np.arange(len(coord1))[match]
    idx2 = idx[match]
    return idx1, idx2, d2d[match], d3d[match]

def _get_join_funcs(colname:str, seplimit:units.Quantity)->dict:
    return {colname:join_skycoord(seplimit, _custom_match_func)}

def merge_photom_tables(tab1:Table, tab2:Table, tol:units.Quantity=1*u.arcsec,
                        table_names:list=None)->Table:
    """
    Given two photometry source catalogs, cross-match and merge them
    using the custom_match_function. This ensures there is a unique
    match between tables as opposed to the default join_skycoord
    behavior which matches multiple objects on the right table to
    a source on the left.
    """
    if table_names is not None:
        assert len(table_names)==2, "Invalid number of table names for two tables."
        assert (type(table_names[0])==str)&(type(table_names[1])==str), "Table names should be strings."


    coord1 = SkyCoord(tab1['ra'], tab1['dec'], unit="deg")
    coord2 = SkyCoord(tab2['ra'], tab2['dec'], unit="deg")
    idx, d2d, _ = coord1.match_to_catalog_sky(coord2)
    match = d2d<tol
    matched_tab2 = tab2[idx][match]
    matched_tab1 = tab1[match]

    # tab1 INTERSECTION tab2
    inner_join = hstack([matched_tab1, matched_tab2],
                        table_names=table_names)
    if table_names is None:
        table_names = ['1', '2']
    tab1_coord_cols = ['ra_'+table_names[0],"dec_"+table_names[0]]
    tab2_coord_cols = ['ra_'+table_names[1],"dec_"+table_names[1]]


    inner_join.remove_columns(tab2_coord_cols)
    inner_join.rename_columns(tab1_coord_cols, ['ra', 'dec'])

    not_matched_tab1 = setdiff(tab1, matched_tab1, keys=['ra', 'dec'])
    not_matched_tab2 = setdiff(tab2, matched_tab2, keys=['ra', 'dec'])

    # (tab1 UNION tab2() - (tab1 INTERSECTION tab2)
    outer_join = join(not_matched_tab1, not_matched_tab2,
                      keys=['ra','dec'], join_type='outer', table_names=table_names)

    #Bring it all together
    return vstack([inner_join, outer_join])