""" Created by Yuxin Dong """
from astropy import units as u
from astropy.coordinates import SkyCoord 
import requests
import json
from collections import OrderedDict


# For searches 

# ID of your Bot:
YOUR_BOT_ID=111031

# name of your Bot:
YOUR_BOT_NAME="Avon"

# API key of your Bot:
api_key="fc32024eaf71cad9e5b3880c833e8b83c676996f"




def parse_coord(ra: str | float, dec: str | float) -> SkyCoord | None:
    """
    Parse an RA/Dec pair into a SkyCoord.

    Args:
        ra (str or float): Right ascension; either in degrees or
            sexagesimal (hh:mm:ss) if ``dec`` is sexagesimal too.
        dec (str or float): Declination; either in degrees or
            sexagesimal (dd:mm:ss) if ``ra`` is sexagesimal too.

    Returns:
        astropy.coordinates.SkyCoord or None: The ICRS coordinate, or None
        (after printing an error) if the input cannot be interpreted.

    """
    if (not (is_number(ra) and is_number(dec)) and
        (':' not in ra and ':' not in dec)):
        error = 'ERROR: cannot interpret: {ra} {dec}'
        print(error.format(ra=ra, dec=dec))
        return(None)

    if (':' in str(ra) and ':' in str(dec)):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        error = 'ERROR: Cannot parse coordinates: {ra} {dec}'
        print(error.format(ra=ra,dec=dec))
        return(None)
    
def is_number(num: str | float) -> bool:
    """
    Check whether a value can be converted to a float.

    Args:
        num (str or float): Value to check.

    Returns:
        bool: True if ``float(num)`` works, False if it raises a ValueError.

    """
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)


### sarching for matching transients using TNS API and a specified radius ###

# function for changing data to json format
def format_to_json(source: str) -> list | dict:
    """
    Extract the reply from the JSON text of a TNS API response.

    Args:
        source (str): JSON text of the response (e.g. ``response.text``).

    Returns:
        list or dict: The ``reply`` entry of the ``data`` of the response.

    """
    # change data to json format and return
    parsed = json.loads(source)   
    result = parsed['data']
    #print(result)
    result = parsed['data']['reply']
    return result



# function for search obj from tutorial
def search(json_list: dict | list) -> "requests.Response | list":
  """
  Search the TNS for objects matching the given search parameters.

  Args:
      json_list (dict or list of tuple): Search parameters of the TNS API
          (e.g. ra, dec, radius, units), which are sent as JSON.

  Returns:
      requests.Response or list: The response of the TNS server. If the request
      fails, a list ``[None, <error message>]`` instead.

  """
  try:
    search_url='https://www.wis-tns.org/api/get/search'
    # url for search obj
    #search_url=url+'/search'
    # headers
    headers={'User-Agent':'tns_marker{"tns_id":'+str(YOUR_BOT_ID)+', "type":"bot",'\
             ' "name":"'+YOUR_BOT_NAME+'"}'}
    # change json_list to json format
    json_file=OrderedDict(json_list)
    # construct a dictionary of api key data and search obj data
    search_data={'api_key':api_key, 'data':json.dumps(json_file)}
    # search obj using request module
    response=requests.post(search_url, headers=headers, data=search_data)
    # return response
    return response
  except Exception as e:
    return [None,'Error message : \n'+str(e)]
