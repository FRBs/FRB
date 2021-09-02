"""
Module for extracting sources and performing photometry
with DECam images
"""

import numpy as np, os, glob

from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.convolution import Gaussian2DKernel
from astropy.wcs import WCS
from astropy.table import Table
import matplotlib.pyplot as plt

import sep
from photutils import Background2D
from photutils.background import MADStdBackgroundRMS, SExtractorBackground
from photutils.segmentation import SourceCatalog, SegmentationImage

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

    return source_cat


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
    source_cat = _source_table(data_sub,segm,bkg_img,1/np.sqrt(wtmap),**sourcecat_kwargs)

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
        # Write
        source_cat.write(os.path.join(outdir, tilename+"_source_cat.dat"), format="ascii.fixed_width", overwrite=True)

        np.savez_compressed(os.path.join(outdir, tilename+"_segmap.npz"), segmap=segmap)
    
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
