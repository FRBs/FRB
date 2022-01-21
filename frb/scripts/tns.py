""" Created by Yuxin Dong """
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from frb.surveys import tns_util

from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Do a TNS Search!')
    parser.add_argument("ra", type=float, help="RA (deg)")
    parser.add_argument("dec", type=float, help="Dec (deg)")
    parser.add_argument("radius", type=float, help="search radius (deg)")
    parser.add_argument("units", type=str, help="units of search radius [arcsec, arcmin, deg]")
    parser.add_argument("--skip", type=str, help="Comma separate list of transient types to skip")
    parser.add_argument("--outfile", type=str, help="Write the results to this file")
    #parser.add_argument("--override", default=False, action='store_true',
    #                    help="Over-ride errors (as possible)? Not recommended")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):

    # coords of your object
    #ra = float(input("Enter RA of your object in deg.: "))
    #dec = float(input("Enter Dec of your object in deg.: "))

    pos = SkyCoord(ra=pargs.ra, dec=pargs.dec, 
                   unit=(u.hourangle,u.deg))

    # SEARCHING

    # ID of your Bot:
    YOUR_BOT_ID=111031

    # name of your Bot:
    YOUR_BOT_NAME="Avon"

    # API key of your Bot:
    api_key="fc32024eaf71cad9e5b3880c833e8b83c676996f"

    # how large of an area to search?
    #radius = float(input("Enter radius to search in: "))
    #units = str(input("Enter the radius unit: ")) #arcsecond, arcmin, deg


    # search obj
    #ra = ra
    #dec = dec
    result = []
    search_obj = [("ra", pargs.ra), ("dec", pargs.dec), 
                  ("radius", pargs.radius), ("units",pargs.units), ]
    response = tns_util.search(search_obj)


    json_data = tns_util.format_to_json(response.text)
    result.append(json_data)
    #print(result)

    new_result = []; new_indices = []
    for i, val in enumerate(result):
        for j, v in enumerate(val):
            if pargs.skip is not None:
                if v['prefix'] in pargs.skip.split(','):
                    continue
            # Grab it
            new_result.append(v)
            new_indices.append(i)

    # write the result out
    df = pd.DataFrame(new_result) 
    if pargs.outfile is not None:
        df.to_csv(pargs.outfile, index=True)

    # Print to screen
    print(df)