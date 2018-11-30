""" I/O related to surveys """
import os


def save_plt(plt, out_dir, root, verbose=False, ftype='png'):
    """
    Save a matplotlib object to disk

    Args:
        plt: matplotlib.pyplot
        out_dir: str
          Folder for output
        root: str
          Root name of the output file
        verbose: bool, optional
        ftype: str
          File type, e.g.  png, pdf

    Returns:

    """
    # Prep
    basename = root+'.{:s}'.format(ftype)
    outfile = os.path.join(out_dir, basename)

    # Write
    plt.savefig(outfile, dpi=300)
    if verbose:
        print("Wrote: {:s}".format(outfile))
        
        
def write_catalog(tbl, out_dir, ftype='ecsv', verbose=False):
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
        verbose: bool, optional

    Returns:

    """
    # Check
    if ftype not in ['ecsv']:
        raise IOError("Unallowed file type: {:s}".format(ftype))
    #
    root = tbl.meta['survey']
    # Outfile
    basename = root+'.{:s}'.format(ftype)
    outfile = os.path.join(out_dir, basename)

    # Write
    tbl.write(outfile, overwrite=True)
    if verbose:
        print("Wrote: {:s}".format(outfile))
