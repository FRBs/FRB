import numpy as np
import matplotlib.pyplot as plt
import math 
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord 
import os
import requests
import json
from collections import OrderedDict
from astropy.coordinates import Angle


# For searches 

# ID of your Bot:
YOUR_BOT_ID=111031

# name of your Bot:
YOUR_BOT_NAME="Avon"

# API key of your Bot:
api_key="fc32024eaf71cad9e5b3880c833e8b83c676996f"




def parse_coord(ra, dec):
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
    
def is_number(num):
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)


### sarching for matching transients using TNS API and a specified radius ###

# function for changing data to json format
def format_to_json(source):
    
    # change data to json format and return
    parsed = json.loads(source)   
    result = parsed['data']
    #print(result)
    result = parsed['data']['reply']
    return result



# function for search obj from tutorial
def search(json_list):
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
