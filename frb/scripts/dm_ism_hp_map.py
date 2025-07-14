import argparse
from frb.dm.dm_ism_models import dm_ism_healpix_map
import healpy as hp

DEFAULT_MAPFILE = dm_ism_healpix_map.get_current_mapfile()


def add_create_map_options(parser=None, usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler='resolve')
    parser.add_argument('--nside', type=int, default=64,
                        help='HEALPix resolution parameter (nside). Default is 64.')
    parser.add_argument('--n_cores', type=int, default=15,
                        help='Number of CPU cores to use for parallel processing. Default is 15.')
    parser.add_argument('--save_path', type=str, default=DEFAULT_MAPFILE,
                        help=f'Path to save the HEALPix map file. Default is "{DEFAULT_MAPFILE}".')
    return parser



def add_get_map_options(parser=None, usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler='resolve')

    parser.add_argument('--mapfile', type=str, default=DEFAULT_MAPFILE,
                        help=f'Path to the HEALPix map file. Default is "{DEFAULT_MAPFILE}".')
    return parser

def add_grab_dm_ism_options(parser=None, usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler='resolve')

    parser.add_argument('--l', type=float, required=True,
                        help='Galactic longitude in degrees.')
    parser.add_argument('--b', type=float, required=True,
                        help='Galactic latitude in degrees.')
    parser.add_argument('--mapfile', type=str, default=DEFAULT_MAPFILE,
                        help=f'Path to the HEALPix map file. Default is "{DEFAULT_MAPFILE}".')

    return parser

def add_plot_map_options(parser=None, usage=None):
    if parser == None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler='resolve')

    parser.add_argument('--mapfile', type=str, default=DEFAULT_MAPFILE,
                        help=f'Path to the HEALPix map file. Default is "{DEFAULT_MAPFILE}".')
    parser.add_argument('--title', type=str, default=r'$DM_{ISM}$ Map',
                        help='Title for the Mollweide plot.')
    parser.add_argument('--min', type=float, default=0,
                        help='Minimum value for the color scale.')
    parser.add_argument('--max', type=float, default=1000,
                        help='Maximum value for the color scale.')

    return parser



def main():
    parser = argparse.ArgumentParser(description='Process DM ISM HEALPix maps.')
    subparsers = parser.add_subparsers(dest='command', required=True)

    # Subcommands
    add_create_map_options(subparsers.add_parser('create', help='Create a HEALPix map.'))
    add_get_map_options(subparsers.add_parser('get', help='Get DM ISM value from the map.'))
    add_grab_dm_ism_options(subparsers.add_parser('grab', help='Grab DM ISM value from the map.'))
    add_plot_map_options(subparsers.add_parser('plot', help='Plot the HEALPix map.'))

    args = parser.parse_args()

    if args.command == 'create':
        dm_ism_healpix_map.create_ne2001_dm_healpix_map(nside=args.nside, n_cores=args.n_cores, save_path=args.save_path)
        print(f"Map saved to {args.save_path}")

    elif args.command == 'get':
        dm_map = dm_ism_healpix_map.get_dm_map(mapfile=args.mapfile)
        nside = hp.get_nside(dm_map)
        print(f"Loaded DM map from {args.mapfile}")
        print(f"Number of pixels: {len(dm_map)}")
        print(f"Nside: {nside}")
        print(f"Min DM: {dm_map.min():.2f}, Max DM: {dm_map.max():.2f}")

    elif args.command == 'grab':
        dm_map = dm_ism_healpix_map.get_dm_map(mapfile=args.mapfile)
        dm_ism = dm_ism_healpix_map.grab_dm_ism_from_healpix_map(args.l, args.b, dm_map)
        print(f"DM ISM value at (l={args.l}, b={args.b}): {dm_ism:.2f} pc cm^-3")

    elif args.command == 'plot':
        dm_ism_healpix_map.plot_mollwiede_view_dm_ism(mapfile=args.mapfile, title=args.title, min=args.min, max=args.max)

    return None


if __name__ == '__main__':
    main()

# This script is designed to be run from the command line.
# It provides functionality to create a HEALPix map of DM ISM values,
# retrieve DM ISM values from the map, grab specific DM ISM values based on coordinates,
# and plot the HEALPix map in a Mollweide projection.
# The script uses argparse for command-line argument parsing.