"""VISTA catalog"""

import re
import warnings
from urllib.parse import urljoin

import numpy as np
from astropy import io, utils
from astropy import units
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils
from frb.galaxies.defs import VISTA_bands

# Dependencies
try:
    from pyvo.dal import sia
except ImportError:
    print("Warning:  You need to install pyvo to retrieve VISTA images")
    _svc = None
else:
    _DEF_ACCESS_URL = "https://datalab.noao.edu/sia/vhs_dr5"
    _svc = sia.SIAService(_DEF_ACCESS_URL)

try:
    import requests
except ImportError:
    requests = None

_VSA_GETIMAGE_FORM_URL = "http://vsa.roe.ac.uk:8080/vdfs/VgetImage_form.jsp"
_VSA_GETIMAGE_ACTION = "./GetImage"
_VSA_ARCHIVE = "VSA"
_VHS_PROGRAMME_ID = "110"
_VISTA_FILTER_IDS = {
    "Z": "1",
    "Y": "2",
    "J": "3",
    "H": "4",
    "KS": "5",
}

# Define the data model for DES data
photom = {}
photom['VISTA'] = {}
photom['VISTA']['VISTA_ID'] = 'sourceid'
photom['VISTA']['ra'] = 'ra2000'
photom['VISTA']['dec'] = 'dec2000'
photom['VISTA']['VISTA_CLASS'] = 'mergedclass' #Class flag,1|0|-1|-2|-3|-9=gal|noise|star|probStar|probGal|saturated
for band in VISTA_bands:
    photom['VISTA']['VISTA_{:s}'.format(band)] = '{:s}petromag'.format(band.lower())
    photom['VISTA']['VISTA_{:s}_err'.format(band)] = '{:s}petromagerr'.format(band.lower())

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['VISTA'] = {'VISTA_ID': int, 'VISTA_CLASS': int}



class VISTA_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the VISTA (VHS) survey

    Child of DL_Survey which uses datalab to access NOAO for catalogs.
    Images are retrieved from the VISTA Science Archive (VSA).

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.dlsurvey.DL_Survey`
            (e.g. ``verbose``)

    """

    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'VISTA'
        self.bands = VISTA_bands
        self.svc = _svc
        self.qc_profile = "default"
        self.database = "vhs_dr5.vhs_cat_v3"

    def _parse_cat_band(self, band: str) -> tuple[list[str], list[str], str]:
        """
        Internal method to generate the bands for grabbing
        a cutout image

        For VISTA, nothing much is necessary.

        Args:
            band (str): Band desired

        Returns:
            tuple: Three items:

                - list of str: Names of the image table columns to select on
                - list of str: Values the columns must have
                - str: Band string for cutout

        """
        table_cols = ['proctype','prodtype']
        col_vals = ['Stack','image']

        return table_cols, col_vals, band

    def _gen_cat_query(self, query_fields: list[str] | None = None,
                       qtype: str = 'main') -> str:
        """
        Generate SQL Query for catalog search

        self.query is modified in place

        Args:
            query_fields (list of str, optional):  Override the default list for the SQL query
            qtype (str, optional):  Type of query to generate.  Currently only 'main' is supported

        Returns:
            str: The SQL query (also stored in ``self.query``)

        Raises:
            IOError: If ``qtype`` is not 'main' and ``query_fields`` is None.

        """
        if query_fields is None:
            query_fields = []
            # Main query
            if qtype == 'main':
                for key,value in photom['VISTA'].items():
                    query_fields += [value]
                database = self.database
            else:
                raise IOError("Bad qtype")

        self.query = dlsurvey._default_query_str(query_fields, database,self.coord,self.radius)

        # Because they HAD to include the epoch in the colname.
        self.query = self.query.replace('ra,dec,','ra2000,dec2000,')
        # Return
        return self.query

    def get_catalog(self, query: str | None = None,
                    query_fields: list[str] | None = None,
                    print_query: bool = False, system: str = 'AB',
                    **kwargs) -> Table:
        """
        Grab a catalog of sources around the input coordinate to the search radius

        Args:
            query (str, optional): SQL query. If None, it is generated
                from ``query_fields``.
            query_fields (list of str, optional): Over-ride list of items to query
            print_query (bool, optional): Print the SQL query generated
            system (str, optional): Magnitude system ['AB', 'Vega']
            **kwargs: Passed to :meth:`frb.surveys.dlsurvey.DL_Survey.get_catalog`
                (e.g. ``timeout``)

        Returns:
            astropy.table.Table:  Catalog of sources returned, with VISTA
            photometry in AB magnitudes by default.  Can be empty.

        Raises:
            RuntimeError: If ``system`` is not 'AB' or 'Vega'.

        """
        # Main DES query
        if query==None:
            self.query = self._gen_cat_query(query_fields=query_fields)
        else:
            self.query = query
        main_cat = super(VISTA_Survey, self).get_catalog(query=self.query, print_query=print_query,
                                                         photomdict=photom['VISTA'],**kwargs)
        if len(main_cat) == 0:
            main_cat = catalog_utils.clean_cat(main_cat, photom['VISTA'], mask_photometry=True)
            main_cat = catalog_utils.ensure_empty_schema(main_cat, list(photom['VISTA'].keys()),
                                                       dtypes=schema_dtypes['VISTA'])
            return main_cat
        # Convert to AB mag
        if system == 'AB':
            #http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=Paranal&gname2=VISTA
            fnu0 = {'VISTA_Y':2087.32,
                    'VISTA_J':1554.03,
                    'VISTA_H':1030.40,
                    'VISTA_Ks':674.83}
            for filt in fnu0.keys():
                main_cat[filt] -= 2.5*np.log10(fnu0[filt]/3630.7805)
        elif system == 'Vega':
            pass
        else:
            raise RuntimeError("Photometry system must be one of 'AB' and 'Vega'")
        
        main_cat = catalog_utils.clean_cat(main_cat, photom['VISTA'], mask_photometry=True)
        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog


    def _select_best_img(self, imgTable: Table, verbose: bool,
                         timeout: int | float = 120) -> io.fits.HDUList:
        """
        Select the best band for a cutout

        Args:
            imgTable (astropy.table.Table): Table of images, with 'exptime'
                and 'access_url' columns
            verbose (bool):  Print status
            timeout (int or float, optional):  How long to wait before timing out, in seconds

        Returns:
            astropy.io.fits.HDUList: The downloaded image with the longest exposure time

        """
        row = imgTable[np.argmax(imgTable['exptime'].data.data.astype('float'))] # pick image with longest exposure time
        url = row['access_url'].decode()
        if verbose:
            print ('downloading deepest stacked image...')

        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat

    @staticmethod
    def _extract_select_map(html: str) -> dict[str, list[str]]:
        """
        Extract select names and their option values from a form page.

        Args:
            html (str): HTML of the form page.

        Returns:
            dict: Maps the name of each select element with options to the
            list of str option values (or labels, if they have no values).

        """
        select_map = {}
        select_pattern = re.compile(r'<select[^>]*name=["\']([^"\']+)["\'][^>]*>(.*?)</select>',
                                    re.IGNORECASE | re.DOTALL)
        option_pattern = re.compile(r'<option[^>]*(?:value=["\']([^"\']*)["\'])?[^>]*>(.*?)</option>',
                                    re.IGNORECASE | re.DOTALL)

        for name, body in select_pattern.findall(html):
            values = []
            for value, text in option_pattern.findall(body):
                token = (value or text or "").strip()
                if token:
                    values.append(token)
            if values:
                select_map[name] = values
        return select_map

    @staticmethod
    def _extract_form_action_method(html: str) -> tuple[str | None, str]:
        """
        Extract form action and method from the first form in page HTML.

        Args:
            html (str): HTML of the form page.

        Returns:
            tuple: ``(action, method)``: the form action (str, None if the
            form has none) and the lower-case method (str, 'post' if the form
            has none).

        """
        action = None
        method = "post"
        form_match = re.search(r'<form[^>]*>', html, flags=re.IGNORECASE)
        if form_match:
            form_tag = form_match.group(0)
            action_match = re.search(r'action=["\']([^"\']+)["\']', form_tag, flags=re.IGNORECASE)
            method_match = re.search(r'method=["\']([^"\']+)["\']', form_tag, flags=re.IGNORECASE)
            if action_match:
                action = action_match.group(1).strip()
            if method_match:
                method = method_match.group(1).strip().lower() or "post"
        return action, method

    @staticmethod
    def _extract_input_names(html: str) -> list[str]:
        """
        Extract all input names from a form page.

        Args:
            html (str): HTML of the form page.

        Returns:
            list of str: The unique names of the input elements, in order.

        """
        input_pattern = re.compile(r'<input[^>]*name=["\']([^"\']+)["\']', re.IGNORECASE)
        return list(dict.fromkeys(input_pattern.findall(html)))

    @staticmethod
    def _extract_select_options(html: str, select_name: str) -> list[tuple[str, str]]:
        """
        Extract option values and labels for a named select element.

        Args:
            html (str): HTML of the form page.
            select_name (str): Name of the select element.

        Returns:
            list of tuple: ``(value, label)`` of str for each option;
            empty if there is no such select element.

        """
        pattern = re.compile(
            rf'<select[^>]*name=["\']{re.escape(select_name)}["\'][^>]*>(.*?)</select>',
            re.IGNORECASE | re.DOTALL,
        )
        match = pattern.search(html)
        if not match:
            return []

        options = []
        for opt in re.finditer(r'<option[^>]*>(.*?)</option>', match.group(1), re.IGNORECASE | re.DOTALL):
            opt_tag = opt.group(0)
            value_match = re.search(r'value=["\']?([^"\'\s>]+)', opt_tag, re.IGNORECASE)
            value = value_match.group(1).strip() if value_match else ""
            label = re.sub(r'<[^>]+>', '', opt.group(1)).strip()
            options.append((value, label))
        return options

    @staticmethod
    def _extract_fits_links(html: str, base_url: str) -> list[str]:
        """
        Extract absolute FITS links from an HTML response.

        Args:
            html (str): HTML of the response.
            base_url (str): URL against which relative links are resolved.

        Returns:
            list of str: The unique absolute URLs of FITS files (or of the
            VSA wrappers of them).

        """
        links = []
        href_pattern = re.compile(r'href=["\']([^"\']+)["\']', re.IGNORECASE)
        text_url_pattern = re.compile(r'https?://[^\s"\'<>]+', re.IGNORECASE)
        fits_ext = (".fits", ".fit", ".fits.fz", ".fit.fz")

        for candidate in href_pattern.findall(html):
            lower = candidate.lower()
            if any(ext in lower for ext in fits_ext) or "getimage.cgi" in lower or "getfimage.cgi" in lower:
                links.append(urljoin(base_url, candidate))

        for candidate in text_url_pattern.findall(html):
            lower = candidate.lower()
            if any(ext in lower for ext in fits_ext) or "getimage.cgi" in lower or "getfimage.cgi" in lower:
                links.append(candidate)

        return list(dict.fromkeys(links))

    def _resolve_vsa_download_link(self, session: "requests.Session", link: str,
                                   timeout: int | float = 120,
                                   verbose: bool = False) -> str:
        """
        Resolve getImage.cgi wrapper links to direct FITS download links.

        Args:
            session (requests.Session): Session used to fetch the wrapper page.
            link (str): URL returned by the VSA.
            timeout (int or float, optional): Seconds to wait for the server.
            verbose (bool, optional): Print status.

        Returns:
            str: The direct download link, or ``link`` itself if it is
            not a wrapper or could not be resolved (a warning is issued
            in the latter case).

        """
        lower = link.lower()
        if "getfimage.cgi" in lower:
            return link
        if "getimage.cgi" not in lower:
            return link

        try:
            wrapper = session.get(link, timeout=timeout)
            wrapper.raise_for_status()
            nested_links = self._extract_fits_links(wrapper.text, wrapper.url)
            if verbose:
                print(f"Resolved wrapper link into {len(nested_links)} nested candidate(s).")
            if not nested_links:
                return link

            for candidate in nested_links:
                if "getfimage.cgi" in candidate.lower():
                    return candidate
            return nested_links[0]
        except Exception as exc:
            warnings.warn(f"Failed to resolve VSA wrapper link: {exc}")
            return link

    @staticmethod
    def _pick_vhs_database(database_options: list[tuple[str, str]]) -> str:
        """
        Pick latest VHS release value from VSA database options.

        Args:
            database_options (list of tuple): ``(value, label)`` of str of the
                database options of the VSA form.

        Returns:
            str: Value of the latest VHS data release, or 'VHSDR7' if there
            are no options.

        """
        if not database_options:
            return "VHSDR7"

        preferred = []
        for value, label in database_options:
            val = (value or "").strip()
            text = (label or "").strip()
            if not val or val.lower() == "none":
                continue
            if val.upper().startswith("VHSDR"):
                suffix = val.upper().replace("VHSDR", "")
                try:
                    rank = int(suffix)
                except ValueError:
                    rank = -1
                preferred.append((rank, val, text))

        if preferred:
            preferred.sort(reverse=True)
            return preferred[0][1]

        for value, _ in database_options:
            val = (value or "").strip()
            if val and val.lower() != "none":
                return val
        return "VHSDR7"

    @staticmethod
    def _to_sexagesimal_strings(coord: SkyCoord) -> tuple[str, str]:
        """
        Convert ICRS coordinates to VSA-friendly sexagesimal strings.

        Args:
            coord (astropy.coordinates.SkyCoord): The coordinate.

        Returns:
            tuple of str: RA (hh:mm:ss.ss) and Dec (+dd:mm:ss.ss).

        """
        ra_str = coord.ra.to_string(unit=units.hourangle, sep=':', precision=2, pad=True)
        dec_str = coord.dec.to_string(unit=units.deg, sep=':', precision=2, pad=True, alwayssign=True)
        return ra_str, dec_str

    @staticmethod
    def _choose_option(values: list[str], contains: str) -> str | None:
        """
        Choose first value containing token, case-insensitive.

        Args:
            values (list of str): Candidate values.
            contains (str): Token to look for.

        Returns:
            str or None: The first value that contains the token, or None
            if there is none.

        """
        token = contains.lower()
        for value in values:
            if token in value.lower():
                return value
        return None

    def _build_vsa_payload(self, html: str, coord: SkyCoord, size_arcmin: float,
                           band: str) -> dict[str, str]:
        """
        Build a permissive form payload from parsed fields and heuristics.

        Args:
            html (str): HTML of the VSA form page.
            coord (astropy.coordinates.SkyCoord): Center of the image.
            size_arcmin (float): Size of the image in arcmin.
            band (str): VISTA band.

        Returns:
            dict: The form fields and their values.

        """
        select_map = self._extract_select_map(html)
        input_names = self._extract_input_names(html)
        payload = {}

        band_lower = band.lower()
        ra_deg = f"{coord.ra.deg:.8f}"
        dec_deg = f"{coord.dec.deg:.8f}"
        size_str = f"{size_arcmin:.6f}"

        for name, values in select_map.items():
            lowered = name.lower()
            selected = values[0]

            if any("j2000" in value.lower() for value in values):
                selected = self._choose_option(values, "j2000") or selected
            elif any(value.lower() == band_lower for value in values):
                selected = band
            elif any("all" == value.lower() for value in values):
                selected = self._choose_option(values, "all") or selected

            if "wave" in lowered or "filter" in lowered or "band" in lowered:
                selected = self._choose_option(values, band_lower) or selected
            elif "coord" in lowered or "system" in lowered:
                selected = self._choose_option(values, "j2000") or selected
            elif "frame" in lowered and any("tilestack" in value.lower() for value in values):
                selected = self._choose_option(values, "tilestack") or selected
            elif "obs" in lowered and any("object" in value.lower() for value in values):
                selected = self._choose_option(values, "object") or selected
            elif ("survey" in lowered or "prog" in lowered) and any("vhs" in value.lower() for value in values):
                selected = self._choose_option(values, "vhs") or selected

            payload[name] = selected

        for name in input_names:
            lowered = name.lower()
            if lowered in payload:
                continue
            if "ra" in lowered and "frame" not in lowered:
                payload[name] = ra_deg
            elif "dec" in lowered:
                payload[name] = dec_deg
            elif ("x" in lowered and "size" in lowered) or lowered in {"xsize", "xs"}:
                payload[name] = size_str
            elif ("y" in lowered and "size" in lowered) or lowered in {"ysize", "ys"}:
                payload[name] = size_str
            elif "multiframe" in lowered or "frameset" in lowered:
                payload[name] = ""
            elif "submit" in lowered:
                payload[name] = "Submit"

        # Extra fallback aliases in case form field names differ from guessed names.
        payload.update({
            "ra": ra_deg,
            "dec": dec_deg,
            "xsize": size_str,
            "ysize": size_str,
            "waveband": band,
            "filter": band,
            "coordSystem": "J2000",
            "frameType": payload.get("frameType", "tilestack"),
            "obsType": payload.get("obsType", "object"),
        })

        return payload

    def _query_vsa_cutout_links(self, imsize: units.Quantity, band: str,
                                timeout: int | float = 120,
                                verbose: bool = False) -> list[str]:
        """
        Query the VSA getImage form and extract candidate FITS links.

        Args:
            imsize (astropy.units.Quantity): Angular size of the image.
            band (str): VISTA band.
            timeout (int or float, optional): Seconds to wait for the server.
            verbose (bool, optional): Print status.

        Returns:
            list of str: URLs of candidate FITS images; empty (with a warning)
            if the query fails or the band is unknown.

        """
        if requests is None:
            warnings.warn("requests is required for VSA image retrieval but is not installed.")
            return []

        size_arcmin = float(imsize.to(units.arcmin).value)
        coord = self.coord.icrs
        ra_str, dec_str = self._to_sexagesimal_strings(coord)
        band_key = band.strip().upper()
        filter_id = _VISTA_FILTER_IDS.get(band_key)
        if filter_id is None:
            warnings.warn(f"No VSA filter mapping found for VISTA band '{band}'.")
            return []

        try:
            with requests.Session() as session:
                # Load form with VHS programme pre-selected so the database list is populated.
                form_params = {
                    "database": "",
                    "programmeID": _VHS_PROGRAMME_ID,
                    "ra": ra_str,
                    "dec": dec_str,
                    "sys": "J",
                    "filterID": filter_id,
                    "xsize": f"{size_arcmin:.6f}",
                    "ysize": f"{size_arcmin:.6f}",
                    "obsType": "object",
                    "frameType": "tilestack",
                    "mfid": "",
                    "fsid": "",
                }
                form_resp = session.get(_VSA_GETIMAGE_FORM_URL, params=form_params, timeout=timeout)
                form_resp.raise_for_status()

                db_opts = self._extract_select_options(form_resp.text, "database")
                database_value = self._pick_vhs_database(db_opts)

                payload = {
                    "archive": _VSA_ARCHIVE,
                    "programmeID": _VHS_PROGRAMME_ID,
                    "database": database_value,
                    "ra": ra_str,
                    "dec": dec_str,
                    "sys": "J",
                    "filterID": filter_id,
                    "xsize": f"{size_arcmin:.6f}",
                    "ysize": f"{size_arcmin:.6f}",
                    "obsType": "object",
                    "frameType": "tilestack",
                    "mfid": "",
                    "fsid": "",
                }

                submit_url = urljoin(form_resp.url, _VSA_GETIMAGE_ACTION)
                response = session.post(submit_url, data=payload, timeout=timeout)

                response.raise_for_status()
                links = self._extract_fits_links(response.text, response.url)
                links = [self._resolve_vsa_download_link(session, link, timeout=timeout, verbose=verbose)
                         for link in links]
                if verbose:
                    print(f"VSA returned {len(links)} cutout link(s) for band {band} ({database_value}).")
                return links
        except Exception as exc:
            warnings.warn(f"VSA query failed for VISTA image retrieval: {exc}")
            return []

    @staticmethod
    def _select_best_vsa_link(links: list[str], band: str) -> str | None:
        """
        Select a deterministic best link, preferring URLs that mention band.

        Args:
            links (list of str): Candidate URLs.
            band (str): VISTA band.

        Returns:
            str or None: The first link that mentions the band, otherwise
            the first link; None if there are no links.

        """
        if not links:
            return None
        band_lower = band.lower()
        preferred = [link for link in links if band_lower in link.lower()]
        return preferred[0] if preferred else links[0]

    def get_image(self, imsize: units.Quantity, band: str | None = None,
                  timeout: int | float = 120, verbose: bool = False
                  ) -> io.fits.PrimaryHDU | None:
        """
        Retrieve a VISTA FITS image through the VSA getImage service.

        Args:
            imsize (astropy.units.Quantity): Angular size of the image.
            band (str, optional): VISTA band (case-insensitive). If None,
                the first band of the survey is used.
            timeout (int or float, optional): Seconds to wait for the server.
            verbose (bool, optional): Print status.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None (with a
            warning) if none could be retrieved.

        Raises:
            TypeError: If ``band`` is not one of the VISTA bands.

        """
        if band is None:
            band = self.bands[0]
            warnings.warn(f"Retrieving VISTA image in default {band} band.")

        allowed = [item.lower() for item in self.bands]
        if band.lower() not in allowed:
            raise TypeError("Allowed filters (case-insensitive) for {:s} photometric bands are {}".format(
                self.survey, self.bands
            ))

        links = self._query_vsa_cutout_links(imsize=imsize, band=band, timeout=timeout, verbose=verbose)
        best_link = self._select_best_vsa_link(links, band)
        if best_link is None:
            warnings.warn(f"No VSA FITS image available for VISTA at requested position in {band} band.")
            return None

        try:
            filename = utils.data.download_file(best_link, cache=True, show_progress=False, timeout=timeout)
            with io.fits.open(filename) as hdul:
                primary = hdul[0]
                if primary.data is not None:
                    return io.fits.PrimaryHDU(data=primary.data, header=primary.header)

                for ext in hdul[1:]:
                    if getattr(ext, "data", None) is not None:
                        return io.fits.PrimaryHDU(data=ext.data, header=ext.header)

                return io.fits.PrimaryHDU(header=primary.header)
        except Exception as exc:
            warnings.warn(f"Failed to download/open VSA FITS image for VISTA: {exc}")
            return None

