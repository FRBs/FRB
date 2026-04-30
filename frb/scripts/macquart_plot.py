"""
Script to generate DM vs z plots for FRBs from the F4 Watchlist. Download all sheets of the watchlist into a folder. Then, run this script with the location of the folder as the --table-loc argument. The script will collate all the tables, calculate DM_ISM, DM_EG and DM_cosmic_est for each FRB, and save the collated table to a CSV file. It will then generate a DM vs z plot (DM_EG or DM_cosmic_est vs z) and save it as a PNG file. You can customize the plot with various command-line arguments. See parse_args() for details.

Typical usage:

python macquart_plot.py --table-loc ../DESI/FRB_tables/ --outfile watchlist_collated_with_script.csv --fig-outfile fig_macquart.png   --color pink --emoji-file Cherry-Blossom-Emoji.png --emoji-zoom 0.02   --plot-dm-eg --plot-running-mean --show-plot
"""

import argparse
import os

import numpy as np

from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from frb.frb import FRB
from frb import mw
from frb.dm.igm import average_DM

from scipy.stats import binned_statistic

from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astropy.cosmology import Planck18


import glob

import re

def parse_jname(coord_str:str)->SkyCoord:
    """
    Use regular expressions to detect instances
    of JXXXXXX.XX+XXXXXX.XX and convert them to
    SkyCoord objects.
    Args:
        coord_str (str): strng containing a
            coordinate in the JXXXXXX+XXXXXX format.
            Essentially, the letter 'J' followed byhttps://github.com/AemulusProject/hmf_emulator
            a number, a sign and another number. The
            length of the numbers can be arbitrarily
            long. It is however assumed to have units
            of hourangle and degrees.
    Returns:
        coord (SkyCoord): A SkyCoord object encoding
            the coordinates.
    """
    
    # Define the regex format
    coord_regex = r"J\d+\.?\d+[+-]?\d+\.?\d+" #J(number)(.?)(number?)(+/-)(number)(.?)(number)
    # Find all instances
    coord = re.findall(coord_regex, coord_str)
    # How many instances found?
    if len(coord) == 1:
        coord = coord[0][1:]
        # Parse delimiter
        if "+" in coord:
            delim = "+"
        else:
            delim = "-"
        # RA and DEC
        try:
            ra, dec = coord.split(delim)
        except:
            import pdb; pdb.set_trace()
        ra = ra[:2]+"h"+ra[2:4]+"m"+ra[4:]+"s"
        dec = delim+dec[:2]+"d"+dec[2:4]+"m"+dec[4:]+"s"
        # Convert to SkyCoord and return
        try:
            return SkyCoord(ra,dec)
        except:
            print(f"Uh-oh. Something went wrong when reading the {coord}")
            import pdb; pdb.set_trace()
    # Curently I only care if the string has one
    # instance of a coordinate.
    else:
        return SkyCoord(np.nan, np.nan, unit="deg")


def calc_dmcosmic_est(DM_FRB, DM_ISM, z, return_dmeg=False):
    """ Calculate the estimated cosmic DM from the FRB DM, ISM DM, and redshift
    Return the error too

    Args:
        DM_FRB (float): The observed DM of the FRB in pc/cm^3
        DM_ISM (float): The estimated DM contribution from the Milky Way ISM in pc/cm^3
        z (float): The redshift of the FRB
        return_dmeg (bool, optional): Whether to return the estimated DM contribution from the host galaxy. Defaults to False.
    Returns:
        tuple: Values, error [float, float] or [np.ndarray, np.ndarray]
    """
    avg_DM_host = 186. # pc/cm^3, rest-frame :: This is the median of the log-normal PDF from James+2022c
    #DM_host_err = 70. # pc/cm^3, rest-frame
    avg_DM_MW_halo = 40. # pc/cm^3, rest-frame
    DM_MW_halo_err = 0.2 # Relative error

    # host_mu = 2.18*np.log(10) # From James+2022c
    # host_sig = 0.48*np.log(10) # From James+2022c
    DM_cosmic_est = DM_FRB - DM_ISM - avg_DM_MW_halo
    if not return_dmeg:
        DM_cosmic_est -= avg_DM_host/(1+z)
    
    # DM_MW
    DM_MW = avg_DM_MW_halo + DM_ISM
    # Error
    #DM_cosmic_err = np.sqrt((DM_host_err/(1+z))**2+(0.2*(DM_MW))**2)
    #DM_cosmic_err = 0.2*(DM_MW)

    # Return
    return DM_cosmic_est#, DM_cosmic_err


def process_watchlist_tables(table_loc:str, outfile:str)->Table:
    """
    Process watchlist tables and save the results to a single file.
    
    Args:
        table_loc (str): The location of the watchlist tables.
        outfile (str): The location to save the processed table.
    
    Returns:
        Table: The processed table containing FRB data and calculated DM values.
    """
    # Ingest data
    # Locate files
    tab_files = glob.glob(os.path.join(table_loc, "F4_Watchlist*.csv"))
    tab_files.sort()
    tabs = []

    # loop through all the spreadsheets
    for file in tab_files:
        # Read csv file
        tab = Table.read(file, format="ascii.csv")
        
        # Make table column names uniform.
        
        try:
            tab = tab[['FRB (TNS)', 'z','DM','FRB Coord', 'Host Coord', 'Survey']]
        except KeyError:
            tab = tab[['FRB (TNS)', 'z','DM_FRB', 'FRB Coord', 'Host Coord', 'Survey']]
            tab.rename_column('DM_FRB', 'DM')
        
        # Make column data types uniform
        tab['DM'] = tab['DM'].astype(str)
        tab['z'] = tab['z'].astype('str')
        tab.rename_column('FRB (TNS)', 'FRB')
        tab.rename_column('FRB Coord', 'coord_str')
        tab.rename_column('Host Coord', 'host_coord_str')
        # Add the cleaned table to the list
        tabs.append(tab)
        
    # Stack the tables
    frb_tab = vstack(tabs).filled(-99.)
    # Remove bad redshift entries and convert to floats
    frb_tab = frb_tab[~np.isin(frb_tab['z'], ['','TBD', None])]
    frb_tab['z'] = frb_tab['z'].astype(float)
    frb_tab = frb_tab[frb_tab['z']>0]

    # Parse the J-names to SkyCoord objects
    frb_coord = [parse_jname(coord) for coord in frb_tab['coord_str']]
    frb_tab['ra'] = [coord.ra.value for coord in frb_coord]
    frb_tab['dec'] = [coord.dec.value for coord in frb_coord]
    host_coord = [parse_jname(coord) for coord in frb_tab['host_coord_str']]
    frb_tab['host_ra'] = [coord.ra.value for coord in host_coord]
    frb_tab['host_dec'] = [coord.dec.value for coord in host_coord]

    # Further clean names and DM entries
    for entry in frb_tab:
        if 'FRB' in entry['FRB']:
            entry['FRB'] = entry['FRB'].replace('FRB', '')
        if '~' in entry['DM']: # I mean, come on! :eye_roll:
            entry['DM'] = entry['DM'].replace('~','')
    frb_tab['DM'] = frb_tab['DM'].astype(float)
    frb_tab = frb_tab[frb_tab['DM']>0]
    frb_tab = frb_tab[['FRB', 'ra', 'dec', 'DM', 'z', 'host_ra', 'host_dec', 'Survey']]
    frb_tab

    # Calculate DM ISM
    frb_tab['DM_ISM']  = [mw.ismDM(coord).value for coord in frb_coord]
    # Subsequently, DM_eg
    frb_tab['DM_EG'] = [calc_dmcosmic_est(entry['DM'], entry['DM_ISM'], entry['z'], return_dmeg=True) for entry in frb_tab]
    # And then, DM_cosmic_est
    frb_tab['DM_cosmic_est'] = [calc_dmcosmic_est(entry['DM'], entry['DM_ISM'], entry['z']) for entry in frb_tab]

    frb_tab.write(outfile, format="ascii.csv", overwrite=True)

    return frb_tab

# Prepare plotting helper functions
def set_fontsize(ax,fsz):
    '''
    Parameters
    ----------
    ax : Matplotlib ax class
    fsz : float
      Font size
    '''
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)

    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles=handles, loc="lower right", fontsize=fsz)

def set_mplrc():
    """
    Update mpl.rcParams
    """
    mpl.rcParams['mathtext.default'] = 'it'
    mpl.rcParams['font.size'] = 24
    mpl.rc('font',family='sans-serif')
    mpl.rcParams['text.latex.preamble'] = r'\boldmath'
    mpl.rc('text', usetex=True)

def fig_macquart(frb_tab,fig=None, ax=None,
                 outfile:str="fig_macquart.png",
                 plot_dm_eg:bool=True,
                 plot_running_mean:bool=True,
                 special_frbs:list=[],
                 plot_kwargs:dict={},
                 color='pink',
                 fontsize=24,
                 emoji_file = None,
                 emoji_zoom = 0.02,
                 special_color="blue",
                 show_plot:bool=True):
    if ax is not None:
        pass
    else:
        fig, ax = plt.subplots(1, 1, **plot_kwargs)

    set_fontsize(ax, fontsize)
    # make the scatter plot

    if len(special_frbs)>0:
        selected_frbs = frb_tab[np.isin(frb_tab['FRB'], special_frbs)]
        plot_special = True
        print(selected_frbs)
    else:
        plot_special = False

    if plot_dm_eg:
        x, y = frb_tab['z'], frb_tab['DM_EG']
        if plot_special:
            x_special, y_special = selected_frbs['z'], selected_frbs['DM_EG']
        ax.scatter(x, y, label=r"$\rm Known~FRBs$", c=color)
        if plot_special:
            ax.scatter(x_special, y_special, c=special_color)
        ylabel = r"$\rm DM_{EG}$"
        dmmacquart, zeval = average_DM(z = 2.3, cumul=True)
        #ax.plot(zeval, dmmacquart.value+186/(1+zeval), color="brown", label=r"$\rm \langle DM_{cosmic} \rangle+186/(1+z)$")
        ax.plot(zeval, dmmacquart.value+120/(1+zeval), color="chocolate",ls="--", label=r"$\rm \langle DM_{cosmic} \rangle+120/(1+z)$")
        #ax.plot(zeval, dmmacquart.value, color="peru",ls=":", label=r"$\rm \langle DM_{cosmic}\rangle$")
        #ax.plot(zeval, 0.6*dmmacquart.value, color="sandybrown",ls="-.", label=r"$\rm 0.6\times\langle DM_{cosmic}\rangle$")
    else:
        x, y = frb_tab['z'], frb_tab['DM_cosmic_est']
        if plot_special:
            x_special, y_special = selected_frbs['z'], selected_frbs['DM_cosmic_est']
        ax.scatter(x, y, label=r"$\rm Known~FRBs$", c=color) #+100/(1+frb_tab['z']),
        if plot_special:
            ax.scatter(x_special, y_special, c=special_color)
        ylabel = r"$\rm DM_{Cosmic}$"
        dmmacquart, zeval = average_DM(z = 2.3, cumul=True)
        ax.plot(zeval, dmmacquart, color="red",ls="--", label=r"$\rm \langle DM_{cosmic} \rangle$")
        #ax.plot(zeval, 0.2*dmmacquart.value, color="sandybrown",ls="-.", label=r"$0.2\times\langle DM_{cosmic}\rangle$")
    
    # Axes labels and limits
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r"$\rm Redshift$")

    ax.set_ylim(0,)

    # add emojis as markers if desired
    if emoji_file is not None:
        emoji_image = plt.imread(emoji_file)
        imagebox = OffsetImage(emoji_image, zoom=emoji_zoom)
        
        for xi, yi in zip(x, y):
            ab = AnnotationBbox(imagebox, (xi, yi), frameon=False, zorder=1)
            ax.add_artist(ab)

    # Plot running mean if desired
    if plot_running_mean:

        # Calculate binned means, stds and counts
        binned_avg, bin_edges, _ = binned_statistic(x, y, statistic='mean', bins=10)
        binned_std, bin_edges, _ = binned_statistic(x, y, statistic='std', bins=10)
        counts, bin_edges, _ =binned_statistic(x, y, statistic='count', bins=10)

        # Calculate bin centers and sizes for plotting
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
        bin_sizes = np.diff(bin_edges)/2

        # Plot
        ax.errorbar(bin_centers, binned_avg,
                    xerr=bin_sizes,
                    yerr=binned_std/np.sqrt(counts), # Error in binned mean = std/sqrt(N)
                    color='r', fmt="o",
                    markersize=5, zorder=99, markeredgecolor="k",
                    label=r"$\rm Binned~mean$")


    plt.legend(loc="best")

    if show_plot:
        plt.show()

    if outfile is not None:
        plt.savefig(outfile, dpi=300)
    
    return fig, ax

def parse_args():
    """Parse command-line arguments for table processing and plotting."""
    parser = argparse.ArgumentParser(
        description="Generate a Macquart-style DM-z plot from F4 watchlist CSV tables."
    )
    parser.add_argument(
        "--table-loc",
        default="./",
        help=(
            "Directory containing watchlist CSVs named like "
            "F4_watchlist*.csv (default: ./)."
        ),
    )
    parser.add_argument(
        "--outfile",
        default="watchlist_collated.csv",
        help="Path to collated watchlist table CSV (default: watchlist_collated.csv).",
    )
    parser.add_argument(
        "--force-rebuild",
        action="store_true",
        help="Rebuild the collated table even if --outfile already exists.",
    )
    parser.add_argument(
        "--fig-outfile",
        default="fig_macquart_with_script.png",
        help="Output filename for the figure (default: fig_macquart_with_script.png).",
    )
    parser.add_argument(
        "--color",
        default="pink",
        help="Scatter color for the FRB points (default: pink).",
    )
    parser.add_argument(
        "--emoji-file",
        default=None,
        help="Path to emoji image for custom point markers (default: frb_emoji.png).",
    )
    parser.add_argument(
        "--emoji-zoom",
        type=float,
        default=0.02,
        help="Zoom level for emoji markers (default: 0.02).",
    )
    parser.add_argument(
        "--font-size",
        type=float,
        default=17,
        help="Final axis/legend font size (default: 17).",
    )
    parser.add_argument(
        "--special-frbs",
        nargs="*",
        default=[],
        help="Optional list of FRB names to highlight.",
    )
    parser.add_argument(
        "--special-color",
        default="blue",
        help="Color used for highlighted special FRBs (default: blue).",
    )
    parser.add_argument(
        "--plot-dm-eg",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Plot DM_EG instead of DM_cosmic_est (default: True).",
    )
    parser.add_argument(
        "--plot-running-mean",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Overplot binned running mean with errors (default: True).",
    )
    parser.add_argument(
        "--show-plot",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Display the figure interactively (default: False).",
    )
    return parser.parse_args()

def main():
    args = parse_args()

    set_mplrc()

    if (not os.path.exists(args.outfile)) or args.force_rebuild:
        frb_tab = process_watchlist_tables(args.table_loc, args.outfile)
    else:
        frb_tab = Table.read(args.outfile, format="ascii.csv")

    _, ax = fig_macquart(
        frb_tab,
        plot_dm_eg=args.plot_dm_eg,
        plot_running_mean=args.plot_running_mean,
        special_frbs=args.special_frbs,
        color=args.color,
        fontsize=args.font_size,
        plot_kwargs={"figsize": (10, 6)},
        emoji_file=args.emoji_file,
        emoji_zoom=args.emoji_zoom,
        special_color=args.special_color,
        show_plot=args.show_plot,
        outfile=args.fig_outfile,
    )


if __name__ == "__main__":
    main()