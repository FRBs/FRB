#!/bin/env python3
import json
import argparse
import urllib.request, urllib.error, urllib.parse
import time
import sys
import csv
import getpass
import os
import os.path
import re
import ssl
from io import StringIO
from . import surveycoord
from . import catalog_utils
from pandas import read_csv
from astropy.table import Table
from .defs import HSC_API_URL as api_url



# adapted from https://hsc-gitlab.mtk.nao.ac.jp/ssp-software/data-access-tools/-/blob/master/pdr3/hscReleaseQuery/hscReleaseQuery.py
version = 20190514.1


# Define the data model for HSC data
photom = {}
photom['HSC'] = {}
HSC_bands = ['g', 'r', 'i', 'z', 'Y']
for band in HSC_bands:
    photom['HSC']['HSC_{:s}'.format(band)] = '{:s}_kronflux_mag'.format(band.lower())
    photom['HSC']['HSC_{:s}_err'.format(band)] = '{:s}_kronflux_magerr'.format(band.lower())
    photom['HSC']['HSC_{:s}_extendedness'.format(band)] = '{:s}_extendedness_value'.format(band.lower())

photom['HSC']['HSC_ID'] = 'object_id'
photom['HSC']['ra'] = 'ra' # r band is the reference band
photom['HSC']['dec'] = 'dec'
photom['HSC']['photo_z'] = 'photoz_best'
photom['HSC']['photo_z_err'] = 'photoz_std_best'

class HSC_Survey(surveycoord.SurveyCoord):
    """
    Class to handle queries on the HSC database

    Args:
        coord (SkyCoord): CoordinAte for surveying around
        radius (Angle): Search radius around the coordinate

    """
    def __init__(self, coord, radius, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        #
        self.survey = 'HSC'
        self.data_release = 'pdr3'


    def get_catalog(self, query_fields=None, query=None, max_time=120,
                    print_query=False, query_table='pdr3_wide.summary',
                    photoz_table = 'mizuki'):
        """
        Query HSC for all objects within a given
        radius of the input coordinates.

        Args:
            query_fields: list, optional
              Column names to be queried. Default values are
              list(photom['HSC'].values()) if None is passed.
            query: str, optional
              Full query as a string to be passed to the database.
              Overrides the default query.
            max_time: float, optional
              The maximum time interval to wait between query status checks. Defaults to 120s.
            print_query: bool, optional
              Print the SQL query for the photo-z values
            query_table: str, optional
              The table to query. Defaults to 'pdr3_wide.forced'

        Returns:
            catalog: astropy.table.Table
              Contains all measurements retieved
              *WARNING* :: The SDSS photometry table frequently has multiple entries for a given
              source, with unique objid values

        """
        if query_fields is None:
            query_fields = list(photom['HSC'].values())
        # Call


        # Now query for photo-z
        if query is None:
            query = f"SELECT {','.join(query_fields)}\n"
            query += f"FROM {query_table}\n"
            iswide = query_table.split(".")[0].split("_")[-1]=='wide'
            if iswide:
                query += f"FULL OUTER JOIN {query_table.split('.')[0]}.photoz_{photoz_table} USING (object_id)"
            query += "WHERE\n"
            query += f"conesearch(coord, {self.coord.ra.value}, {self.coord.dec.value}, {self.radius.to('arcsec').value})"

        if print_query:
            print(query)

        # SQL command
        query_cat = run_query(query, max_time=max_time,
                              release_version=self.data_release, delete_job=True)

        catalog = catalog_utils.clean_cat(query_cat, photom['HSC'])

        self.catalog = catalog_utils.sort_by_separation(catalog, self.coord, radec=('ra','dec'), add_sep=True)

        # Meta
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        # Validate
        self.validate_catalog()

        # Return
        return self.catalog.copy()
        
class QueryError(Exception):
    pass
    
def run_query(query:str,
              user:str=None,
              release_version:str='pdr3',
              preview:bool=False,
              out_format:str='csv',
              delete_job:bool=False,
              max_time:int=120
              ):
    """
    Submits a query to the HSC database and downloads the results in the specified format.

    Args:
        query (str): The SQL query to submit to the HSC database.
        user (str, optional): The account name to use for authentication. Defaults to None.
        release (str, optional): The release version of the HSC database to query. Defaults to 'pdr3'.
        preview (bool, optional): Whether to use quick mode (short timeout). Defaults to False.
        out_format (str, optional): The format in which to download the query results. Defaults to 'csv'.
        delete_job (bool, optional): Whether to delete the job after downloading the results. Defaults to False.
        max_time (int, optional): The maximum time interval to wait for checking query status. Defaults to 120s.

    Raises:
        urllib.error.HTTPError: If there is an HTTP error while submitting the query.
        QueryError: If there is an error with the query itself.

    Returns:
        None
    """
    user, password = getCredentials()
    credential = {'account_name': user, 'password': password}
    sql = query

    job = None

    try:
        if preview:
            preview(credential, sql, sys.stdout)
        else:
            job = submitJob(credential, sql,
                            out_format=out_format,
                            release_version=release_version)
            blockUntilJobFinishes(credential, job['id'],
                                  max_time=max_time)
            res = download(credential, job['id'])
            pseudo_file = StringIO(res.read().decode('utf-8').split("# ")[1])
            table = Table.from_pandas(read_csv(pseudo_file)).filled(-99.)
            if delete_job:
                deleteJob(credential, job['id'])
            return table
    except urllib.error.HTTPError as e:
        if e.code == 401:
            print('invalid id or password.', file=sys.stderr)
        if e.code == 406:
            print(e.read(), file=sys.stderr)
        else:
            print(e, file=sys.stderr)
    except QueryError as e:
        print(e, file=sys.stderr)
    except KeyboardInterrupt:
        if job is not None:
            jobCancel(credential, job['id'])
        raise


def httpJsonPost(url, data):
    data['clientVersion'] = version
    postData = json.dumps(data)
    return httpPost(url, postData, {'Content-type': 'application/json'})


def httpPost(url, postData, headers):
    req = urllib.request.Request(url, postData.encode('utf-8'), headers)
    res = urllib.request.urlopen(req)
    return res


def submitJob(credential, sql,
              out_format:str="csv", release_version:str="pdr3"):
    url = api_url + 'submit'
    catalog_job = {
        'sql'                     : sql,
        'out_format'              : out_format,
        'include_metainfo_to_body': False,
        'release_version'         : release_version,
    }
    postData = {'credential': credential, 'catalog_job': catalog_job, 'nomail': True, 'skip_syntax_check': False}
    res = httpJsonPost(url, postData)
    job = json.load(res)
    return job


def jobStatus(credential, job_id):
    url = api_url + 'status'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    job = json.load(res)
    return job


def jobCancel(credential, job_id):
    url = api_url + 'cancel'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


def preview(credential, sql, out, release_version:str="pdr3"):
    url = api_url + 'preview'
    catalog_job = {
        'sql'             : sql,
        'release_version' : release_version,
    }
    postData = {'credential': credential, 'catalog_job': catalog_job}
    res = httpJsonPost(url, postData)
    result = json.load(res)

    writer = csv.writer(out)
    # writer.writerow(result['result']['fields'])
    for row in result['result']['rows']:
        writer.writerow(row)

    if result['result']['count'] > len(result['result']['rows']):
        raise QueryError('only top %d records are displayed !' % len(result['result']['rows']))


def blockUntilJobFinishes(credential, job_id, max_time=120):
    interval = 1
    while True:
        time.sleep(interval)
        job = jobStatus(credential, job_id)
        if job['status'] == 'error':
            raise QueryError('query error: ' + job['error'])
        if job['status'] == 'done':
            break
        interval *= 2
        if interval > max_time:
            interval = max_time


def download(credential, job_id):
    url = api_url + 'download'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    return res


def deleteJob(credential, job_id):
    url = api_url + 'delete'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


def getCredentials():
    password_from_envvar = os.environ.get("HSC_SSP_CAS_PASSWORD")
    user_from_envvar = os.environ.get("HSC_SSP_CAS_USER")
    if isinstance(user_from_envvar, str) & isinstance(password_from_envvar, str):
        return user_from_envvar, password_from_envvar
    else:
        raise ValueError("Please set the environment variables HSC_SSP_CAS_USER and HSC_SSP_CAS_PASSWORD to your CAS credentials. Follow the instructions at https://hsc-release.mtk.nao.ac.jp/doc/index.php/data-access__pdr3/ to register.")

