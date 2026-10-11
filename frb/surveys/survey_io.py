""" I/O related to surveys """
import os
from types import ModuleType


def save_plt(plt: ModuleType, out_dir: str, root: str, verbose: bool = False,
             ftype: str = 'png') -> None:
    """
    Save a matplotlib object to disk

    Args:
        plt (module): ``matplotlib.pyplot`` holding the figure to save
        out_dir (str): Folder for output
        root (str): Root name of the output file
        verbose (bool, optional): Print the name of the output file
        ftype (str, optional): File type, e.g.  png, pdf

    Returns:
        None

    """
    # Prep
    basename = root+'.{:s}'.format(ftype)
    outfile = os.path.join(out_dir, basename)

    # Write
    plt.savefig(outfile, dpi=300)
    if verbose:
        print("Wrote: {:s}".format(outfile))
        
        

