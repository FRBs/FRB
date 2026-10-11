""" Parent class for surveying at a given coordinate"""
from abc import ABCMeta

import os
from typing import Any

from astropy import units
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table

from frb.surveys import images
from frb.surveys import survey_io
from frb.surveys import catalog_utils


class SurveyCoord(object):
    """
    Parent class of surveying around an input coordinate.

    API semantics for survey children:
    - ``get_catalog`` returns an astropy Table around ``coord`` within ``radius``.
    - ``get_image`` returns FITS-like image products when survey provides them.
    - ``get_cutout`` returns rendered products (e.g. PNG/JPEG) when available.

    Compatibility policy:
    Surveys that only provide FITS image retrieval should expose that via
    ``get_image``. If those surveys keep ``get_cutout`` for backward
    compatibility, it should call ``get_image`` and emit ``DeprecationWarning``.


    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        verbose (bool, optional): Print status messages by default in
            methods that accept a ``verbose`` argument.

    Attributes:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        verbose (bool): Default verbosity
        survey (str): Name of the survey; set by child classes
        catalog (astropy.table.Table): Catalog of sources; None until
            ``get_catalog`` is run
        cutout (numpy.ndarray): Cutout image; None until one is retrieved
        cutout_size (astropy.units.Quantity): Angular size of the cutout;
            None until a cutout is retrieved

    """

    __metaclass__ = ABCMeta

    def __init__(self, coord: SkyCoord, radius: Angle, verbose: bool = False):
        # Load up
        self.coord = coord
        self.radius = radius
        self.verbose = verbose

        # Typically set items
        self.survey = None

        # Standard products
        self.catalog = None
        self.cutout = None
        self.cutout_size = None

    def get_catalog(self) -> Table:
        """
        Run survey catalog query.

        Child classes should set and return ``self.catalog`` as an astropy Table.
        Expected output contract for normalized survey catalogs:
        - Coordinate columns: ``ra``, ``dec``.
        - Metadata keys: ``radius``, ``survey``.
        - Separation column: ``separation`` in arcmin.

        Returns:
            astropy.table.Table: Survey catalog.

        """
        pass

    def get_cutout(self, imsize: units.Quantity) -> Any:
        """
        Retrieve rendered cutout product (e.g. PNG/JPEG).

        For FITS products, use ``get_image`` as canonical method. FITS-only
        surveys may keep ``get_cutout`` as deprecated alias to ``get_image``
        for backward compatibility.

        Args:
            imsize (astropy.units.Quantity): Angular size of desired cutout.

        Returns:
            object or None: Rendered image-like product, depending on survey.
            The base class returns None.

        """
        return None

    def get_image(self, imsize: units.Quantity, band: str | None = None) -> Any:
        """
        Retrieve FITS-like image product.

        Args:
            imsize (astropy.units.Quantity): Angular size of desired image.
            band (str, optional): Filter/band identifier if required by survey.

        Returns:
            object: FITS HDU or survey-specific FITS-like image product.
            The base class returns None.

        """
        pass

    def validate_catalog(self) -> None:
        """
        Validate and normalize catalog to enforce uniform output contract.
        
        Ensures: lowercase ra/dec, required metadata (radius, survey), 
        and separation column in arcmin for all survey catalogs.
        ``self.catalog`` is replaced by the normalized catalog.
        Does nothing if ``self.catalog`` is None.

        Returns:
            None

        Raises:
            AssertionError: If the normalized catalog lacks the 'ra' or 'dec'
                columns or the 'radius' or 'survey' metadata.

        """
        if self.catalog is None:
            return
        
        # Normalize catalog to enforce uniform contract
        self.catalog = catalog_utils.normalize_catalog(
            self.catalog, 
            self.coord, 
            self.radius, 
            self.survey,
            add_sep=True
        )
        
        # Validate minimum required structure
        assert 'ra' in self.catalog.keys(), "Normalized catalog missing 'ra' column"
        assert 'dec' in self.catalog.keys(), "Normalized catalog missing 'dec' column"
        assert 'radius' in self.catalog.meta.keys(), "Catalog missing 'radius' metadata"
        assert 'survey' in self.catalog.meta.keys(), "Catalog missing 'survey' metadata"
            
    def write_catalog(self, out_dir: str, ftype: str = 'ecsv',
                      verbose: bool | None = None, create_dirs: bool = False,
                      overwrite: bool = True) -> None:
        """
        Write ``self.catalog`` to disk as ``<out_dir>/<self.survey>.<ftype>``

        Args:
            out_dir (str): Folder for output
            ftype (str, optional): File type. Only 'ecsv' is allowed.
            verbose (bool, optional): Print the name of the output file.
                If None, uses ``self.verbose``.
            create_dirs (bool, optional): Create the folders to the output dir
                (if needed)?
            overwrite (bool, optional): Overwrite the existing file?

        Returns:
            None

        Raises:
            IOError: If ``ftype`` is not allowed.

        """
        if verbose is None:
            verbose = self.verbose
        # Check
        if ftype not in ['ecsv']:
            raise IOError("Unallowed file type: {:s}".format(ftype))
        #
        root = self.survey

        # Generate output folder?
        if create_dirs:
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

        # Outfile
        basename = root+'.{:s}'.format(ftype)
        outfile = os.path.join(out_dir, basename)

        if (not overwrite) and (os.path.isfile(outfile)):
            print("Output catalog file already exists.  Use overwrite=True as desired")
            return

        # Write
        self.catalog.write(outfile, overwrite=overwrite)
        if verbose:
            print("Wrote: {:s}".format(outfile))

    def write_cutout(self, output_dir: str = './', root: str | None = None,
                     verbose: bool | None = None) -> None:
        """
        Write the cutout image to disk

        Args:
            output_dir (str, optional): Folder for output
            root (str, optional): Root name of the output file.
                Defaults to ``<self.survey>_cutout``.
            verbose (bool, optional): Print the name of the output file.
                If None, uses ``self.verbose``.

        Returns:
            None

        """
        if root is None:
            root = self.survey+'_cutout'
        if verbose is None:
            verbose = self.verbose
        if self.cutout is None:
            print("Need to get the cutout image first!  Use get_cutout()")
        # Prep plot
        plt = images.gen_snapshot_plt(self.cutout, self.cutout_size)
        survey_io.save_plt(plt, output_dir, root, verbose=verbose)
