""" Simple utilities for figures"""

import numpy as np
import matplotlib as mpl


def log_me(val, err):
    """
    Generate log and error from linear input
    
    Args:
        val (float): 
        err (float): 

    Returns:
        float, (float/None):
            Returns none if the err is negative

    """
    if err < 0.:
        xerr = None
    else:
        xerr = np.array([[np.log10(val) - np.log10(val - err)],
                     [-np.log10(val) + np.log10(val + err)]])
    return np.log10(val), xerr


def set_fontsize(ax,fsz):
    """
    Set the fontsize throughout an Axis
    
    Args:
        ax (Matplotlib Axis): 
        fsz (float): Font size

    Returns:

    """
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)


def set_bokeh_fontsize(p, fsz):
    """
    Adjust font size for Bokeh axes

    Args:
        p (Bokeh plot class):
        sz (int): Font size
    """
    p.xaxis.axis_label_text_font_size = '{:d}pt'.format(fsz)
    p.xaxis.major_label_text_font_size = "{:d}pt".format(fsz)
    #
    p.yaxis.axis_label_text_font_size = '{:d}pt'.format(fsz)
    p.yaxis.major_label_text_font_size = "{:d}pt".format(fsz)


def set_mplrc():
    """
    Font fussing for matplotlib

    Returns:

    """
    mpl.rcParams['mathtext.default'] = 'it'
    mpl.rcParams['font.size'] = 12
    mpl.rc('font',family='Times New Roman')
    #mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
    mpl.rc('text', usetex=True)
