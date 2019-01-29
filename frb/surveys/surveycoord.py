""" Parent class for surveying at a given coordinate"""
from abc import ABCMeta

import os

from frb.surveys import images
from frb.surveys import survey_io


class SurveyCoord(object):
    """
    Parent class of surveying around an input coordinate

    See the children for specific methods

    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate

    """

    __metaclass__ = ABCMeta

    def __init__(self, coord, radius, verbose=False):
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

    def get_catalog(self):
        """

        Returns:
            self.catalog

        """
        pass

    def get_cutout(self, imsize):
        return None

    def get_image(self, imsize, filter):
        pass

    def validate_catalog(self):
        if len(self.catalog) > 0:
            # Columns
            assert 'ra' in self.catalog.keys()
            assert 'dec' in self.catalog.keys()
            # Meta
            assert 'radius' in self.catalog.meta.keys()
            assert 'survey' in self.catalog.meta.keys()
            
    def write_catalog(self, out_dir, ftype='ecsv', verbose=None, create_dirs=False,
                      overwrite=True):
        """
        Write an input astropy Table to disk

        Args:
            tbl: astropy.table.Table
            out_dir: str
              Folder for output
            root: str
              Root name of the output file
            ftype: str, optional
              File type, e.g. ecsv
            create_dirs: bool, optional
              Create the folders to the output dir (if needed)?
            overwrite: bool, optional
              Overwrite the existing file?
            verbose: bool, optional

        Returns:

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

    def write_cutout(self, output_dir='./', root=None, verbose=None):
        """
        Write the cutout image to disk

        Args:
            output_dir: str
            root: str, optional
            verbose: bool, optional

        Returns:

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
