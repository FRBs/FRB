""" A script to crossmatch TNS transients with KKO FRBs and plot the 2D Gaussian elliptical distribution
@author: Yuxin Dong
last edited: Feb 12, 2025 """
import numpy as np
import pandas as pd
import json
import time
import requests
import os
import sys
from collections import OrderedDict

from astropy.coordinates import SkyCoord
import ligo.skymap.plot # KEEP: needed for projections in matplotlib
from astropy import units as u

from scipy.stats import chi2
from astropy.io import ascii

from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import matplotlib


import argparse

api_key = os.environ['TNS_API_KEY']
tns = "sandbox.wis-tns.org" # sandbox
url_tns_api = "https://" + tns + "/api"

def parser():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Script to run TNS queries on FRBs.')

    parser.add_argument("--single", action='store_true', help="Indicate if the query is for a single FRB object.")

    parser.add_argument("--filename", type=str, help="a list of FRBs for a batch query.")
    parser.add_argument("--name", type=str, help="FRB TNS name of the single FRB object.")
    parser.add_argument("--ra", type=float, help="Right Ascension of the single FRB object.")
    parser.add_argument("--dec", type=float, help="Declination of the single FRB object.")
    parser.add_argument("--theta", type=float, help="Theta value (make sure you add a negative sign) of the single FRB object.")
    parser.add_argument("--a", type=float, help="Semi-major axis for the single FRB object.")
    parser.add_argument("--b", type=float, help="Semi-minor axis for the single FRB object.")

    parser.add_argument("--radius", type=float, default=3.0, help="Radius for the query (default is 3.0) in arcmin.")

    return parser.parse_args()


def check_tns_api_keywords():
    for key in ['TNS_BOT_ID','TNS_BOT_NAME','TNS_API_KEY']:
        if key not in os.environ.keys():
            raise Exception(f'Add {key} to your environmental variables.')


 # this is for getting metadata
def set_bot_tns_marker():
    bot_id = os.environ['TNS_BOT_ID']
    bot_name = os.environ['TNS_BOT_NAME']
    tns_marker = 'tns_marker{"tns_id": "' + bot_id + '", "type": "bot", "name": "' + bot_name + '"}'
    return tns_marker
    
# New
def search(json_list, url_tns_api):
    try:
        search_url = url_tns_api + "/get/search" #'https://www.wis-tns.org/api/get/search'
        tns_marker = set_bot_tns_marker()
        headers = {'User-Agent': tns_marker}
        json_file = OrderedDict(json_list)
        search_data = {'api_key': api_key, 'data': json.dumps(json_file)}

        #print("Sending request to:", search_url)
        #print("Headers:", headers)
        #print("Request JSON:", search_data)

        response = requests.post(search_url, headers=headers, data=search_data)

        #print("Status Code:", response.status_code)
        #print("Response Text:", response.text)
        return response  # Ensure this always returns a Response object
    
    except Exception as e:
        print("Error occurred:", str(e))
        return None  # Return None instead of a list


def format_to_json(source):
    #parsed = json.loads(source, object_pairs_hook = OrderedDict)
    #result = json.dumps(parsed, indent = 4) this reverse it back to a string instead of json, so don't do that
    return json.loads(source, object_pairs_hook = OrderedDict)


def tns_query(ra, dec, radius, frb_name, units='arcmin', initial_delay=10, max_delay=500, outfile='query_results_final.txt'):

    '''
   Queries transients in TNS at the FRB position with a specfied radius.

    Parameters: 
    -----------
    ra (float): right ascension of the FRB in deg
    dec (float): declination of the FRB in deg
    radius (float): search radius in arcmin
    frb_name (str): TNS name of the FRB
    
    Returns: 
    --------
    query output: dict
        a dictionary of transients found within the search radius near an FRB position
    '''

    search_obj = [("ra", ra), ("dec", dec), ("radius", radius), ("units", units)]
    max_retries = 10
    delay = initial_delay
    attempt = 0
    results_dict = {}
    while True:
        response = search(search_obj, url_tns_api)
        if response.status_code == 429:
            if attempt < max_retries:
                print(f"Throttled. Retrying in {delay} seconds... (Attempt {attempt + 1}/{max_retries})")
                time.sleep(delay)
                delay = min(delay * 2, max_delay)  # Exponential backoff
                attempt += 1
            else:
                print("All retry attempts failed. Unable to get a successful response.")
                break
        else:
            response.raise_for_status()
            result = format_to_json(response.text)
            print(result)
            if result:
                data_list = result.get('data', []) # we just want the data portion of the result
                for item in data_list:
                    if item['prefix'] == 'FRB':
                            print(f"Skipping object with FRB prefix: {item['objname']}")
                            continue 
                    obj_id = item['objid']
                    results_dict[obj_id] = {
                        'FRB Name': frb_name,
                        'Object Name': item['objname'],
                        'Prefix': item['prefix'],
                        'Object ID': obj_id
                    }
                print(f"Results for {frb_name} have been added to the dictionary.")
            else: 
                print(f"No Results found for {frb_name}")
            break
            # if you want to write them into a file
            #if result:
            #    with open(outfile, 'a') as f:
            #        for item in result:
            #            f.write(f"FRB Name: {frb_name}, Object Name: {item['objname']}, Prefix: {item['prefix']}, Object ID: {item['objid']}\n")
            #    print(f"Reply content for {frb_name} has been written to {outfile}")
            #    return outfile
    return results_dict



# The new get_metadata
def get(objname, url_tns_api):
    get_url = url_tns_api + "/get/object"
    tns_marker = set_bot_tns_marker()
    headers = {'User-Agent': tns_marker}
    get_obj = [("objname", objname), ("objid", "")]
    headers = {'User-Agent': tns_marker}
    json_file = OrderedDict(get_obj)
    get_data = {'api_key': api_key, 'data': json.dumps(json_file)}
    response = requests.post(get_url, headers = headers, data = get_data)
    if response.status_code == 200:
        return json.loads(response.text, object_pairs_hook=OrderedDict)
    else:
        print(f"Error fetching metadata for {objname}: {response.status_code}")
    return response


def read_final_catalog(filename):
    f = ascii.read(filename)

    #filtered_f = f[f['a_err'] > 0.00056] 
    #filtered_f = f[(f['a_err'] > 0.00055556) & (f['b_err'] > 0.00055556)] # larger than 2 arcseconds

    
    name = np.array(f['name'].data)
    ra = np.array(f['ra_frb'].data)
    dec = np.array(f['dec_frb'].data)
    theta = -np.array(f['theta'].data) #need to flip that sign
    b = np.array(f['b_err'].data)
    a = np.array(f['a_err'].data)
    #DM = filtered_f['DM'].data
  
    return name, ra, dec, theta, a, b


def cov_matrix(a, b, theta):

    """
    Calculate the covariance matrix for a 2D ellipse.
    
    Parameters:
    a (float): Semi-major axis of the ellipse
    b (float): Semi-minor axis of the ellipse
    theta (float): Position angle of the ellipse in degrees
    
    Returns:
    np.ndarray: 2x2 covariance matrix
    """

    # Convert theta to radians
    theta_rad = np.radians(theta + 90)  # Adjust angle for correct orientation
    # Rotation matrix for an ellipse with position angle
    R = np.array([[np.cos(theta_rad), -np.sin(theta_rad)], 
                  [np.sin(theta_rad), np.cos(theta_rad)]])
    # Diagonal matrix with square of semi-major and semi-minor axes
    D = np.diag([a**2, b**2])
    # Covariance matrix
    covariance_matrix = R @ D @ R.T
    return covariance_matrix


def mahalanobis_distance(point, frbcenter, cov_matrix):
    """
    Calculate the Mahalanobis distance between a given point and the center.
    effectively the Z-score in 1D, and it can be a proxy for how many sigma away you are from the mean.

    Parameters:
    point (array-like): The point for which to calculate the distance (should be [RA, Dec]).
    frbcenter (array-like): The center point, usually the mean [RA, Dec].
    cov_matrix (array-like): The covariance matrix of the data.

    Returns:
    float: The Mahalanobis distance.
    """

    # the frbcenter in this case is the 'mean'
    diff = np.array(point) - np.array(frbcenter)
    
    # correction term for spherical geometry
    correction = diff[0] * (1 / np.cos(np.radians(point[1])))
    diff_corrected = np.array([correction, diff[1]])
    inv_cov_matrix = np.linalg.inv(cov_matrix)
    md = np.sqrt(diff_corrected.T @ inv_cov_matrix @ diff_corrected)
    
    return md


def percentile(mahalanobis_distance, df=2):
    p_value = chi2.cdf(mahalanobis_distance**2, df)
    return p_value

def gauss_contour(frbcenter, cov_matrix, semi_major, transient_name,
                  transient_position=None, levels=[0.68, 0.95, 0.99]):
    
    """
    Plot Gaussian contours around a given FRB center based on its covariance matrix.
    
    This function computes and visualizes the Gaussian contours (ellipses) 
    that represent the uncertainty in the position of the FRB as defined by 
    its covariance matrix. Then, it plots the position of the transient along with 
    the Mahalanobis distance from the FRB center.
    
    Parameters:
    frbcenter (Astropy SkyCoord): The central position of the FRB.
    cov_matrix (array-like): The covariance matrix representing the uncertainties in the FRB's position.
    semi_major (float): The semi-major axis of the Gaussian ellipse in degrees.
    transient_name (str): The name of the transient source.
    transient_position (Astropy SkyCoord): Transient position.
    levels: List of confidence levels for the contours (default: [0.68, 0.95, 0.99]).
    """
    
    eigvals, eigvecs = np.linalg.eigh(cov_matrix)
    order = eigvals.argsort()[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    
    width, height = 2 * np.sqrt(eigvals)
    theta = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
    threshold = chi2.ppf(0.9999, df=2)
    
    # 2 arcmins size, quite arbitrary so can be changed
    buffer = 2 * u.arcsec
    if transient_position is not None:
        size = (frbcenter.separation(transient_position)+buffer).to(u.arcmin)
    else:
        size = 10*u.arcmin #semi_major*60*30

    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(
    projection='astro zoom',
    center=frbcenter,
    radius=size)

    # Plot FRB center
    ax.plot(frbcenter.ra.deg, frbcenter.dec.deg, 
        transform=ax.get_transform('world'), color='black', marker='x', 
        markersize=5)
    
    cmap = matplotlib.colormaps.get_cmap('viridis')

    for i, level in enumerate(levels):
        chi_square_val = chi2.ppf(level, 2)
        color = cmap(i / len(levels))
        ellipse = Ellipse(
            xy=(frbcenter.ra.deg, frbcenter.dec.deg),
            width=width * np.sqrt(chi_square_val),
            height=height * np.sqrt(chi_square_val),
            angle=theta,
            edgecolor=color,
            fc='None',
            lw=2,
            label=f'{int(level*100)}%',
            transform=ax.get_transform('world')
        )
        ax.add_patch(ellipse)
    
    if transient_position is not None:
        ax.plot(transient_position.ra.deg, transient_position.dec.deg, 
                transform=ax.get_transform('world'), 
                color='#fb5607', marker='o', label=transient_name)
        
        
        md = mahalanobis_distance([transient_position.ra.deg, transient_position.dec.deg],
                                  [frbcenter.ra.deg, frbcenter.dec.deg], cov_matrix)
        pt = percentile(md)
        ax.annotate(f'{(1-pt)*100:.3f}%', xy=[transient_position.ra.deg, transient_position.dec.deg],
                    xytext=(12, 12), textcoords='offset points')     
        
    ax.set_xlabel('RA (J2000)', fontsize=14)
    ax.set_ylabel('DEC (J2000)', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(False)
    plt.tight_layout()
    #plt.savefig(f'{transient_name}_gaussian_map.png')


def main(filename, name, ra, dec, theta, a, b, radius, single_obj=False):

    """
    Main function to query the TNS for transient data; can be batched or a single object search.
    Parameters:
    filename: The path to the ASCII file containing the catalog data (only used if single_obj is False).
    name (str): FRB name (only used if single_obj is True).
    ra (float): FRB RA in degrees 
    dec (float): FRB Dec in degrees
    theta: Position Angle in degrees (only used if single_obj is True).
    a: semi-minor axis in degrees (only used if single_obj is True).
    b: semi-minor axis in degrees (only used if single_obj is True).
    radius: The search radius (arcmin).
    single_obj: Boolean indicating whether to perform a query for a single object (True) or multiple (False).

    Returns:
    trans_results: dictionary 
    trans_metadata: dictionary
    2D Gaussian map: png
    """ 
     
    trans_results = {}
    if not single_obj:
        name, ra, dec, theta, a, b = read_final_catalog(filename)
        for n, r, d in zip(name, ra, dec):
            frb_results = tns_query(r, d, frb_name=n, radius=radius)
            trans_results.update(frb_results)
    else:
        frb_results = tns_query(ra, dec, frb_name=name, radius=radius)
        trans_results.update(frb_results)

    #grab all the metadata from TNS for the transients
    trans_metadata = {}
    for obj_id, data in trans_results.items():
        objname = data['Object Name']
        FRBname = data['FRB Name']
        print(f'Fetching metadata for object: {objname}')
        metadata = get(objname, url_tns_api) # grab metadata
        # fill free to add anything other info you want from TNS
        if metadata:
            data_dict = metadata.get('data', {})
            formatted_data = {
                'frbname': FRBname,
                'objname': data_dict.get('objname', ''),
                'prefix': data_dict.get('name_prefix'),
                'type': data_dict.get('object_type', ''),
                'radeg': data_dict.get('radeg', ''),
                'decdeg': data_dict.get('decdeg', ''),
                'ra': data_dict.get('ra', ''),
                'dec': data_dict.get('dec', ''),
                'redshift': data_dict.get('redshift', ''),
                'hostname': data_dict.get('hostname', ''),
               'host_redshift': data_dict.get('host_redshift'),
            }
            trans_metadata[obj_id] = formatted_data

    if single_obj:
        for obj_id, data in trans_metadata.items():
            cov = cov_matrix(a, b, theta)
            transient_pos = SkyCoord(ra=data['radeg'], dec=data['decdeg'], unit='deg')
            frbcenter = SkyCoord(ra=ra, dec=dec, unit='deg')
            gauss_contour(frbcenter, cov, a, data['objname'], transient_pos)
            plt.savefig(f'{name}_{data["objname"]}_gaussian_map.png')
            plt.clf() # Clear the figure for the next plot 
            break  # Exit the loop after finding a match
    else:
        # now we plot the 2D Gaussian maps 
        for n, r, d, t, a, b in zip(name, ra, dec, theta, a, b):
            for obj_id, data in trans_metadata.items():
                if data['frbname'] == n:
                    # semi-minor goes first
                    cov = cov_matrix(a,b,t)
                    transient_pos = SkyCoord(ra=data['radeg'],dec=data['decdeg'],unit='deg')
                    frbcenter = SkyCoord(ra=r,dec=d,unit='deg')
                    gauss_contour(frbcenter, cov, a, data['objname'], transient_pos)
                        
                    # Save the plot as a PNG file
                    plt.savefig(f'{n}_{data["objname"]}_gaussian_map.png')
                    plt.clf()

if __name__ == "__main__":
    # Verify that you've added these to your env var 
    args = parser()
    check_tns_api_keywords()
    if args.single:
        if args.ra is None or args.dec is None or args.theta is None or args.a is None or args.b is None:
            print("Error: Missing parameters for single object query.")
            exit(1)
        matched_transient_data = main(filename=None,
                                      name=args.name, 
                                      ra=args.ra, 
                                      dec=args.dec,
                                      theta=args.theta, 
                                      a=args.a, b=args.b, radius=args.radius, single_obj=True) # in degrees
    else:
        frb_file = args.filename #sys.argv[1]
        matched_transient_data = main(frb_file, name=None, ra=None, dec=None, theta=None, a=None, b=None,
                                      radius=args.radius, single_obj=False)