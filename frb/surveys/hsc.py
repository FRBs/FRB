#!/bin/env python3
import json
import urllib.request, urllib.error
import time
import sys
import csv
import os
import warnings
from io import StringIO
from http.client import HTTPResponse
from typing import TextIO
from astropy.coordinates import Angle, SkyCoord
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
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """
    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        #
        self.survey = 'HSC'
        self.data_release = 'pdr3'


    def get_catalog(self, query_fields: list[str] | None = None,
                    query: str | None = None, timeout: int | float = 120,
                    print_query: bool = False,
                    query_table: str = 'pdr3_wide.summary',
                    photoz_table: str = 'mizuki', **kwargs) -> Table:
        """
        Query HSC for all objects within a given
        radius of the input coordinates.

        Args:
            query_fields (list of str, optional): Column names to be
                queried. Default values are
                list(photom['HSC'].values()) if None is passed.
            query (str, optional): Full query as a string to be passed to
                the database. Overrides the default query.
            timeout (int or float, optional): The maximum time interval to
                wait between query status checks. Defaults to 120s.
            print_query (bool, optional): Print the SQL query
            query_table (str, optional): The table to query. Defaults to
                'pdr3_wide.summary'
            photoz_table (str, optional): Photo-z table joined to the
                query table if it belongs to the 'wide' release.
            **kwargs: Only the deprecated ``max_time`` (use ``timeout``) is
                accepted.

        Returns:
            astropy.table.Table: Contains all measurements retrieved

        Raises:
            TypeError: If an unexpected keyword argument is given, or both
                ``timeout`` and the deprecated ``max_time`` are specified.

        """
        if 'max_time' in kwargs:
            warnings.warn(
                "'max_time' is deprecated; use 'timeout' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if timeout != 120:
                raise TypeError("Specify only one of 'timeout' or deprecated 'max_time'.")
            timeout = kwargs.pop('max_time')
        if kwargs:
            raise TypeError(f"Unexpected keyword arguments: {list(kwargs.keys())}")

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
        query_cat = run_query(query, timeout=timeout,
                              release_version=self.data_release, delete_job=True)

        catalog = catalog_utils.clean_cat(query_cat, photom['HSC'], mask_photometry=True)

        self.catalog = catalog_utils.sort_by_separation(catalog, self.coord, radec=('ra','dec'), add_sep=True)

        # Meta
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        # Validate
        self.validate_catalog()

        # Return
        return self.catalog.copy()
        
class QueryError(Exception):
    """Raised when there is an error in a query to the HSC database
    (or the HSC credentials are missing)."""
    pass
    
def run_query(query: str,
              user: str | None = None,
              release_version: str = 'pdr3',
              preview: bool = False,
              out_format: str = 'csv',
              delete_job: bool = False,
              timeout: int | float = 120,
              max_time: int | float | None = None
              ) -> Table | None:
    """
    Submits a query to the HSC database and downloads the results in the specified format.


    Args:
        query (str): The SQL query to submit to the HSC database.
        user (str, optional): Not used; the credentials are read from the
            environment by :func:`getCredentials`. Defaults to None.
        release_version (str, optional): The release version of the HSC database to query. Defaults to 'pdr3'.
        preview (bool, optional): Whether to use quick mode (short timeout). Defaults to False.
        out_format (str, optional): The format in which to download the query results. Defaults to 'csv'.
        delete_job (bool, optional): Whether to delete the job after downloading the results. Defaults to False.
        timeout (int or float, optional): The maximum time interval to wait for checking query status. Defaults to 120s.
        max_time (int or float, optional): Deprecated; use ``timeout``.

    Raises:
        QueryError: If the HSC credentials are not set in the environment.
        TypeError: If both ``timeout`` and the deprecated ``max_time`` are specified.

    Returns:
        astropy.table.Table or None: The query results, with missing values filled
        with -99. None if the query fails with an HTTP or query error (the
        error is printed to stderr) or ``preview`` is True.
    """
    user, password = getCredentials()
    credential = {'account_name': user, 'password': password}
    sql = query

    if max_time is not None:
        warnings.warn(
            "'max_time' is deprecated; use 'timeout' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if timeout != 120:
            raise TypeError("Specify only one of 'timeout' or deprecated 'max_time'.")
        timeout = max_time

    job = None

    try:
        if preview:
            preview(credential, sql, sys.stdout)
        else:
            job = submitJob(credential, sql,
                            out_format=out_format,
                            release_version=release_version)
            blockUntilJobFinishes(credential, job['id'],
                                  timeout=timeout)
            res = download(credential, job['id'])
            pseudo_file = StringIO(res.read().decode('utf-8').split("# ")[1])
            table = catalog_utils.fill_masked(Table.from_pandas(read_csv(pseudo_file)), -99.)
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


def httpJsonPost(url: str, data: dict) -> HTTPResponse:
    """
    POST a JSON payload to the HSC API.

    Args:
        url (str): URL of the API endpoint.
        data (dict): Payload; the client version is added to it in place.

    Returns:
        http.client.HTTPResponse: The response of the server.

    """
    data['clientVersion'] = version
    postData = json.dumps(data)
    return httpPost(url, postData, {'Content-type': 'application/json'})


def httpPost(url: str, postData: str, headers: dict) -> HTTPResponse:
    """
    POST data to the HSC API.

    Args:
        url (str): URL of the API endpoint.
        postData (str): Payload, which is UTF-8 encoded before sending.
        headers (dict): HTTP headers of the request.

    Returns:
        http.client.HTTPResponse: The response of the server.

    """
    req = urllib.request.Request(url, postData.encode('utf-8'), headers)
    res = urllib.request.urlopen(req)
    return res


def submitJob(credential: dict, sql: str,
              out_format: str = "csv", release_version: str = "pdr3") -> dict:
    """
    Submit a query job to the HSC database.

    Args:
        credential (dict): ``account_name`` and ``password`` for the HSC database.
        sql (str): The SQL query.
        out_format (str, optional): Format of the query results.
        release_version (str, optional): Release of the HSC database to query.

    Returns:
        dict: The job description returned by the server, including its 'id'.

    """
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


def jobStatus(credential: dict, job_id: str) -> dict:
    """
    Get the status of a job on the HSC database.

    Args:
        credential (dict): ``account_name`` and ``password`` for the HSC database.
        job_id (str): ID of the job, as returned by :func:`submitJob`.

    Returns:
        dict: The job description returned by the server, including its 'status'.

    """
    url = api_url + 'status'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    job = json.load(res)
    return job


def jobCancel(credential: dict, job_id: str) -> None:
    """
    Cancel a job on the HSC database.

    Args:
        credential (dict): ``account_name`` and ``password`` for the HSC database.
        job_id (str): ID of the job, as returned by :func:`submitJob`.

    """
    url = api_url + 'cancel'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


def preview(credential: dict, sql: str, out: TextIO, release_version: str = "pdr3") -> None:
    """
    Run a query in preview mode and write the returned rows as CSV.

    Args:
        credential (dict): ``account_name`` and ``password`` for the HSC database.
        sql (str): The SQL query.
        out (file-like): Writable text stream the CSV rows are written to.
        release_version (str, optional): Release of the HSC database to query.

    Raises:
        QueryError: If the preview holds only the top rows of the results.

    """
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


def blockUntilJobFinishes(credential: dict, job_id: str,
                          timeout: int | float = 120,
                          max_time: int | float | None = None) -> None:
    """
    Wait until a job on the HSC database is done.

    Args:
        credential (dict): ``account_name`` and ``password`` for the HSC database.
        job_id (str): ID of the job, as returned by :func:`submitJob`.
        timeout (int or float, optional): The maximum time interval, in seconds,
            to wait between status checks.
        max_time (int or float, optional): Deprecated; use ``timeout``.

    Raises:
        QueryError: If the job ends in an error.
        TypeError: If both ``timeout`` and the deprecated ``max_time`` are specified.

    """
    if max_time is not None:
        warnings.warn(
            "'max_time' is deprecated; use 'timeout' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if timeout != 120:
            raise TypeError("Specify only one of 'timeout' or deprecated 'max_time'.")
        timeout = max_time

    interval = 1
    while True:
        time.sleep(interval)
        job = jobStatus(credential, job_id)
        if job['status'] == 'error':
            raise QueryError('query error: ' + job['error'])
        if job['status'] == 'done':
            break
        interval *= 2
        if interval > timeout:
            interval = timeout


def download(credential: dict, job_id: str) -> HTTPResponse:
    """
    Download the results of a job on the HSC database.

    Args:
        credential (dict): ``account_name`` and ``password`` for the HSC database.
        job_id (str): ID of the job, as returned by :func:`submitJob`.

    Returns:
        http.client.HTTPResponse: The response of the server, holding the results.

    """
    url = api_url + 'download'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    return res


def deleteJob(credential: dict, job_id: str) -> None:
    """
    Delete a job on the HSC database.

    Args:
        credential (dict): ``account_name`` and ``password`` for the HSC database.
        job_id (str): ID of the job, as returned by :func:`submitJob`.

    """
    url = api_url + 'delete'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


def getCredentials() -> tuple[str, str]:
    """
    Read the HSC credentials from the environment.

    Uses the environment variables HSC_SSP_CAS_USER and HSC_SSP_CAS_PASSWORD.

    Returns:
        tuple of str: The user name and the password.

    Raises:
        QueryError: If either environment variable is not set.

    """
    password_from_envvar = os.environ.get("HSC_SSP_CAS_PASSWORD")
    user_from_envvar = os.environ.get("HSC_SSP_CAS_USER")
    if isinstance(user_from_envvar, str) & isinstance(password_from_envvar, str):
        return user_from_envvar, password_from_envvar
    else:
        raise QueryError("Please set the environment variables HSC_SSP_CAS_USER and HSC_SSP_CAS_PASSWORD to your CAS credentials. Follow the instructions at https://hsc-release.mtk.nao.ac.jp/doc/index.php/data-access__pdr3/ to register.")

