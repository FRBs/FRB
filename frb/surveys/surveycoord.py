""" Parent class for surveying at a given coordinate"""
from abc import ABCMeta

from frb.surveys import images
from frb.surveys import survey_io

class SurveyCoord(object):
    """

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
        pass

    def get_cutout(self, imsize):
        pass

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
            
    def write_catalog(self, out_dir='./', verbose=False):
        """
        Write catalog to disk
        Simple wrapper to survey_io.write_catalog

        Args:
            out_dir: str
              Folder for output
            verbose: bool, optional

        Returns:

        """
        survey_io.write_catalog(self.catalog, out_dir, verbose=verbose)

    def write_cutout(self, output_dir='./', root=None, verbose=None):
        if root is None:
            root = self.survey+'_cutout'
        if verbose is None:
            verbose = self.verbose
        if self.cutout is None:
            print("Need to get the cutout image first!  Use get_cutout()")
        # Prep plot
        plt = images.gen_snapshot_plt(self.cutout, self.cutout_size)
        survey_io.save_plt(plt, output_dir, root, verbose=verbose)
