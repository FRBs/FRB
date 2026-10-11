""" Module for image routines"""

from io import BytesIO
from types import ModuleType
from astropy import units
try:
    from PIL import Image
except ImportError:
    print("Warning: You need to install PIL to write SDSS cutout images")

try:
    import requests
except ImportError:
    print("Warning: You need to install requests to handle SDSS images")

from matplotlib import pyplot as plt

def grab_from_url(url: str) -> "Image.Image":
    """
    Grab a PIL Image from a URL

    Args:
        url (str): URL

    Returns:
        PIL.Image.Image: Image retrieved from the URL

    """
    # Simple calls
    rtv = requests.get(url)
    img = Image.open(BytesIO(rtv.content))
    # Return
    return img


def gen_snapshot_plt(img: "Image.Image", imsize: units.Quantity,
                     show: bool = False) -> ModuleType:
    """
    Generate a simple figure from an input PIL.Image

    Args:
        img (PIL.Image.Image): Image to plot
        imsize (astropy.coordinates.Angle or astropy.units.Quantity):
            Angular dimension of the image
        show (bool, optional): Show to the screen?  
           If done, will need to regenerate to then save to disk

    Returns:
        module: ``matplotlib.pyplot``, which allows one to further modify the plot

    """
    # Convert to arcsec and float
    i_arcsec = imsize.to('arcsec').value
    #
    plt.clf()
    plt.imshow(img, aspect='equal', extent=(-i_arcsec / 2., i_arcsec / 2, -i_arcsec / 2., i_arcsec / 2))
    # Label me
    plt.xlabel('Relative arcsec', fontsize=20)
    xpos = 0.22 * i_arcsec
    ypos = 0.02 * i_arcsec
    plt.text(-i_arcsec / 2. - xpos, 0., 'EAST', rotation=90., fontsize=20)
    plt.text(0., i_arcsec / 2. + ypos, 'NORTH', fontsize=20, horizontalalignment='center')
    # Show?
    if show:
        plt.show()
    # Return
    return plt

