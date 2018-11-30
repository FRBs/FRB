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
        
        

