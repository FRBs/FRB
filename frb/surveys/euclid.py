"""
Methods related to Euclid survey queries via astroquery.esa.euclid

Provides access to Euclid Q1 data release through the ESA Euclid Archive,
including photometry (VIS+NIR), morphology, and optional spectroscopy.
"""

import warnings
import signal
import os
import shutil
from collections.abc import Iterator
from contextlib import contextmanager
from types import FrameType
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.table import Table

try:
    from astroquery.esa.euclid import Euclid
    Euclid.ROW_LIMIT = -1  # Remove row limit for queries
except ImportError:
    print("Warning: You need to install astroquery to use the Euclid survey tools...")

from frb.surveys import surveycoord
from frb.surveys import catalog_utils


@contextmanager
def _query_timeout(timeout_seconds: int | float | None) -> Iterator[None]:
    """
    Context manager to enforce a hard timeout on blocking archive queries.

    Uses SIGALRM, so it only works in the main thread on Unix.

    Args:
        timeout_seconds (int or float or None): Seconds after which a
            ``TimeoutError`` is raised. No timeout is enforced if it is None
            or not positive.

    Raises:
        TimeoutError: If the body of the ``with`` block takes longer than
            ``timeout_seconds``.

    """
    if timeout_seconds is None or timeout_seconds <= 0:
        yield
        return

    def _handler(signum: int, frame: FrameType | None) -> None:
        """
        SIGALRM handler that raises a ``TimeoutError``.

        Args:
            signum (int): Signal number.
            frame (types.FrameType or None): Current stack frame.

        """
        raise TimeoutError(f"Euclid query timed out after {timeout_seconds} seconds")

    previous_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _handler)
    signal.setitimer(signal.ITIMER_REAL, float(timeout_seconds))
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)

# Define the data model for Euclid photometry
# Main source: catalogue.mer_catalogue (MER = Multi-band Extraction Reduction)
photom = {}
photom['Euclid'] = {}
Euclid_bands = ['VIS', 'J', 'H', 'Y']
# VIS is detection band, J/H/Y are NISP IR bands
# Map FRB standard filter names to Euclid database columns
photom['Euclid']['Euclid_VIS'] = 'flux_vis_sersic'
photom['Euclid']['Euclid_VIS_err'] = 'fluxerr_vis_sersic'
photom['Euclid']['Euclid_J'] = 'flux_j_sersic'
photom['Euclid']['Euclid_J_err'] = 'fluxerr_j_sersic'
photom['Euclid']['Euclid_H'] = 'flux_h_sersic'
photom['Euclid']['Euclid_H_err'] = 'fluxerr_h_sersic'
photom['Euclid']['Euclid_Y'] = 'flux_y_sersic'
photom['Euclid']['Euclid_Y_err'] = 'fluxerr_y_sersic'
photom['Euclid']['Euclid_ID'] = 'object_id'
photom['Euclid']['ra'] = 'right_ascension'
photom['Euclid']['dec'] = 'declination'

# Optional morphological parameters
photom['Euclid']['Euclid_ellipticity'] = 'ellipticity'
photom['Euclid']['Euclid_kron_radius'] = 'kron_radius'
photom['Euclid']['Euclid_segmentation_area'] = 'segmentation_area'

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['Euclid'] = {'Euclid_ID': int, 'Euclid_segmentation_area': int}

_EUCLID_FLUX_SCALE = 1e-6 / 3630.7805


def _euclid_flux_to_abmag(flux_microjy: np.ndarray,
                          fluxerr_microjy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert Euclid microJy fluxes and errors to AB magnitudes.

    Args:
        flux_microjy (array-like): Fluxes in microJy.
        fluxerr_microjy (array-like): Flux errors in microJy.

    Returns:
        tuple of numpy.ndarray: The AB magnitudes and their errors.
        Both are NaN where the flux is not finite and positive; the error
        is also NaN where the flux error is not finite and positive.

    """
    flux = np.asarray(flux_microjy, dtype=float)
    fluxerr = np.asarray(fluxerr_microjy, dtype=float)
    mag = np.full(flux.shape, np.nan, dtype=float)
    mag_err = np.full(flux.shape, np.nan, dtype=float)

    good_flux = np.isfinite(flux) & (flux > 0)
    mag[good_flux] = -2.5 * np.log10(flux[good_flux] * _EUCLID_FLUX_SCALE)

    good_err = good_flux & np.isfinite(fluxerr) & (fluxerr > 0)
    mag_err[good_err] = (2.5 / np.log(10.0)) * fluxerr[good_err] / flux[good_err]
    return mag, mag_err

# Define the data model for Euclid spectroscopy (if available)
spectrom = {}
spectrom['Euclid'] = {}
spectrom['Euclid']['Euclid_spec_ID'] = 'source_id'
spectrom['Euclid']['Euclid_spec_z'] = 'redshift'
spectrom['Euclid']['ra'] = 'right_ascension'
spectrom['Euclid']['dec'] = 'declination'


class Euclid_Survey(surveycoord.SurveyCoord):
    """
    Class to handle queries on the Euclid survey database.
    
    Queries the ESA Euclid Archive for photometry and optional spectroscopy
    from the Q1 (March 2025) public data release. Primary data source is the
    mer_catalogue table which contains VIS+NIR photometry and morphology.

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    Examples:
        >>> from astropy.coordinates import SkyCoord
        >>> from astropy import units as u
        >>> from frb.surveys.euclid import Euclid_Survey
        >>> coord = SkyCoord(ra=267.78*u.deg, dec=65.53*u.deg)
        >>> radius = 10*u.arcmin
        >>> survey = Euclid_Survey(coord, radius)
        >>> catalog = survey.get_catalog()
    """

    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        self.survey = 'Euclid'

    def get_catalog(self, query_fields: list[str] | None = None,
                    timeout: int | float = 120, check_spectra: bool = False) -> Table:
        """
        Query Euclid for all objects within a given radius of input coordinates.

        Queries the mer_catalogue table for photometry and morphology.
        Optionally searches for associated spectroscopy if photometry is found.

        Args:
            query_fields (list of str, optional):
                List of column names to retrieve from database.
                If None, uses default set: object_id, right_ascension, declination,
                plus all photometric and morphological columns.
            timeout (int or float, optional):
                Query timeout in seconds. Default: 120 s.
            check_spectra (bool, optional):
                If True, check for associated spectra for photometric sources
                (see :meth:`spectra_exist`) and add the boolean column
                'Euclid_has_spectrum'. Default: False.

        Returns:
            astropy.table.Table:
                Catalog of sources with standardized FRB column names.
                Includes ra, dec, and Euclid_{VIS,J,H,Y} magnitudes.
                Empty table if no sources found.

        Raises:
            TimeoutError: If the archive does not respond within ``timeout``.

        """

        # Default fields to query
        if query_fields is None:
            query_fields = [
                'object_id', 'right_ascension', 'declination',
                'flux_vis_sersic', 'fluxerr_vis_sersic',
                'flux_j_sersic', 'fluxerr_j_sersic',
                'flux_h_sersic', 'fluxerr_h_sersic',
                'flux_y_sersic', 'fluxerr_y_sersic',
                'ellipticity', 'kron_radius', 'segmentation_area',
                'vis_det', 'det_quality_flag'
            ]

        # Use cone_search for efficient regional query
        # Note: cone_search returns results with 'dist' column showing angular distance
        if self.verbose:
            print(f"Querying Euclid mer_catalogue around {self.coord.ra.deg:.4f}, {self.coord.dec.deg:.4f}")
            print(f"Search radius: {self.radius.to(u.arcmin).value:.2f} arcmin")

        # Perform cone search
        # By default, Euclid.cone_search uses table_name='catalogue.mer_catalogue'
        with _query_timeout(timeout):
            job = Euclid.cone_search(
                coordinate=self.coord,
                radius=self.radius,
                columns=query_fields,
                async_job=True
            )

        # Get results
        with _query_timeout(timeout):
            photom_catalog = job.get_results()

        if photom_catalog is None or len(photom_catalog) == 0:
            self.catalog = catalog_utils.ensure_empty_schema(
                Table(), list(photom['Euclid'].keys()), dtypes=schema_dtypes['Euclid']
            )
            self.catalog.meta['radius'] = self.radius
            self.catalog.meta['survey'] = self.survey
            if self.verbose:
                print("No sources found in Euclid.")
            self.validate_catalog()
            return self.catalog

        for band in Euclid_bands:
            flux_key = photom['Euclid'][f'Euclid_{band}']
            err_key = photom['Euclid'][f'Euclid_{band}_err']
            if flux_key in photom_catalog.colnames and err_key in photom_catalog.colnames:
                photom_catalog[flux_key], photom_catalog[err_key] = _euclid_flux_to_abmag(
                    photom_catalog[flux_key], photom_catalog[err_key]
                )

        # Clean up catalog - rename columns to FRB standard names
        self.catalog = catalog_utils.clean_cat(photom_catalog, photom['Euclid'], mask_photometry=True)

        # Add metadata
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        if self.verbose:
            print(f"Found {len(self.catalog)} sources in Euclid")

        # Optionally check for associated spectra
        # This is slow so only do if requested and if we have photometric sources to check
        if check_spectra:
            if self.verbose:
                print("Checking for associated spectra...")
            euclid_ids = self.catalog['Euclid_ID']
            has_spec = self.spectra_exist(euclid_ids)
            self.catalog['Euclid_has_spectrum'] = has_spec

        self.validate_catalog()
        return self.catalog.copy()

    def spectra_exist(self, euclid_ids: list | np.ndarray | int) -> np.ndarray:
        """
        Check if spectroscopic data exists for photometric sources.

        This method queries Euclid datalinks using the Euclid object IDs.

        Args:
            euclid_ids (list or numpy.ndarray or int): Euclid object ID(s)
                to check.

        Returns:
            numpy.ndarray: Boolean array, True where there is
            a spectrum for the corresponding ID in ``euclid_ids``.

        """
        euclid_ids = np.atleast_1d(euclid_ids)
        ncat = len(euclid_ids)

        if ncat>100:
            warnings.warn(f"Checking spectra for {ncat} sources may take a long time. Consider limiting to a smaller subset.")

        # Run a get_datalink query
        result = Euclid.get_datalinks(ids = euclid_ids)
        

        if result is not None and len(result) > 0:
            unique_ids = np.unique(result['ID'])
            unique_ids = [uid.replace('sedm ', '') for uid in unique_ids]  # Clean up IDs
        
        has_spec = np.isin(euclid_ids, unique_ids)
        return has_spec

    def get_spectrum(self, euclid_id: int, output_folder: str | None = None,
                     timeout: int | float = 120) -> tuple[None, None] | None:
        """
        Retrieve the spectrum for a given Euclid object ID.

        Args:
            euclid_id (int): The Euclid object ID to retrieve the spectrum for.
            output_folder (str, optional): Output folder for the spectrum FITS
                (written as ``spectrum.fits``). If None, a folder named after
                ``euclid_id`` is used.
            timeout (int or float, optional): Query timeout in seconds. Default: 120 s.

        Returns:
            tuple or None: ``(None, None)`` if no spectrum was found or the
            retrieval failed. If a spectrum was found it is written to disk
            and nothing (None) is returned.

        """

        try:
            if self.verbose:
                print(f"Attempting to retrieve spectrum for Euclid ID {euclid_id}")

            with _query_timeout(timeout):
                has_spec = self.spectra_exist([euclid_id])[0]

            if not has_spec:
                if self.verbose:
                    print(f"No datalinks found for Euclid ID {euclid_id}")
                return None, None
            else:
                if self.verbose:
                    print(f"Datalinks found for Euclid ID {euclid_id}.")

                # Download spectra
                if output_folder is None:
                    output_folder = str(euclid_id) # Use ID as folder name to avoid collisions
                with _query_timeout(timeout):
                    
                    Euclid.get_spectrum(source_id=euclid_id, output_file=f'{output_folder}/spectrum.fits', verbose=False)
        except Exception as e:
            if self.verbose:
                print(f"Spectrum retrieval failed for Euclid ID {euclid_id}: {e}")
            return None, None
                
    def get_image(self, imsize: u.Quantity | None = None,
                  output_file: str | None = None, verbose: bool | None = None,
                  timeout: int | float = 120
                  ) -> tuple[np.ndarray, fits.Header] | tuple[None, None]:
        """
        Get a FITS image cutout of a Euclid MER background-subtracted mosaic image.

        Queries the mosaic_product table to find MER background-subtracted
        mosaics covering the target region, then retrieves a cutout.
        The cutout and its header are also stored in ``self.cutout`` and
        ``self.cutout_hdr``.

        Args:
            imsize (astropy.units.Quantity, optional):
                Size of cutout image. Default: 2 arcmin.
            output_file (str, optional):
                Output filename for cutout FITS. If None, uses temporary location.
            verbose (bool, optional):
                Verbosity. If None, uses self.verbose.
            timeout (int or float, optional):
                Query timeout in seconds. Default: 120 s.

        Returns:
            tuple:
                (data_array, fits_header) of types (numpy.ndarray,
                astropy.io.fits.Header) if cutout successful,
                (None, None) otherwise.

        Note:
            MER cutouts are background-subtracted VIS-band mosaic images.
            Each ~1'x1' cutout weighs ~5.5 MB and takes <1 second to download.
        """
        if verbose is None:
            verbose = self.verbose

        if imsize is None:
            imsize = 2 * u.arcmin


        if verbose:
            print(f"Attempting to retrieve Euclid cutout: {imsize.to(u.arcmin).value:.2f} arcmin")

        # Query for a VIS mosaic that intersects the requested region.
        # Euclid.get_cutout works on MER mosaics, so we fetch one matching file path.
        half_size = 0.5 * imsize
        radius_deg = half_size.to(u.deg).value
        query = (
            "SELECT TOP 1 file_path, file_name, tile_index "
            "FROM q1.mosaic_product "
            "WHERE instrument_name='VIS' "
            f"AND INTERSECTS(CIRCLE({self.coord.ra.deg:.8f}, {self.coord.dec.deg:.8f}, {radius_deg:.8f}), fov)=1 "
            "ORDER BY creation_date DESC"
        )

        if verbose:
            print(query)

        with _query_timeout(timeout):
            mosaic_job = Euclid.launch_job_async(query, verbose=False)
        with _query_timeout(timeout):
            mosaic_table = mosaic_job.get_results()

        if mosaic_table is None or len(mosaic_table) == 0:
            if verbose:
                print("No Euclid VIS mosaic found for requested region")
            self.cutout = None
            self.cutout_hdr = None
            return None, None

        file_path = f"{mosaic_table['file_path'][0]}/{mosaic_table['file_name'][0]}"
        tile_index = str(mosaic_table['tile_index'][0])
        with _query_timeout(timeout):
            cutout_paths = Euclid.get_cutout(
                file_path=file_path,
                instrument='VIS',
                id=tile_index,
                coordinate=self.coord,
                radius=half_size,
                output_file=output_file,
                verbose=False
            )

        if cutout_paths is None or len(cutout_paths) == 0:
            if verbose:
                print("Euclid.get_cutout returned no files")
            self.cutout = None
            self.cutout_hdr = None
            return None, None

        cutout_file = cutout_paths[0]
        with fits.open(cutout_file) as hdul:
            self.cutout = hdul[0].data
            self.cutout_hdr = hdul[0].header.copy()

        # Astroquery writes to a temporary path when output_file is not set.
        # We keep data in memory but remove that on-disk temp artifact.
        if output_file is None:
            try:
                if os.path.isfile(cutout_file):
                    os.remove(cutout_file)
                parent_dir = os.path.dirname(cutout_file)
                if parent_dir and os.path.basename(parent_dir).startswith('temp_'):
                    shutil.rmtree(parent_dir, ignore_errors=True)
            except OSError:
                pass

        self.cutout_size = imsize
            
        return self.cutout, self.cutout_hdr

    def get_cutout(self, imsize: u.Quantity | None = None,
                   output_file: str | None = None, verbose: bool | None = None,
                   timeout: int | float = 120
                   ) -> tuple[np.ndarray, fits.Header] | tuple[None, None]:
        """
        Deprecated alias for get_image().

        Args:
            imsize (astropy.units.Quantity, optional): Size of cutout image.
            output_file (str, optional): Output filename for cutout FITS.
            verbose (bool, optional): Verbosity. If None, uses self.verbose.
            timeout (int or float, optional): Query timeout in seconds.

        Returns:
            tuple: See :meth:`get_image`.

        """
        warnings.warn(
            "get_cutout() returns FITS products for this survey and is deprecated; "
            "use get_image() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.get_image(imsize=imsize, output_file=output_file, verbose=verbose, timeout=timeout)
