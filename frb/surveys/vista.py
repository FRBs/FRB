"""VISTA catalog"""

import re
import warnings
from urllib.parse import urljoin

import numpy as np
from astropy import io, utils
from astropy import units

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



class VISTA_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DECaL survey

    Child of DL_Survey which uses datalab to access NOAO


    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate

    """

    def __init__(self, coord, radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'VISTA'
        self.bands = VISTA_bands
        self.svc = _svc
        self.qc_profile = "default"
        self.database = "vhs_dr5.vhs_cat_v3"

    def _parse_cat_band(self,band):
        """
        Internal method to generate the bands for grabbing
        a cutout image

        For DES, nothing much is necessary.


        Args:
            band (str): Band desired


        Returns:
            list, list, str:  Table columns, Column values, band string for cutout

        """
        table_cols = ['proctype','prodtype']
        col_vals = ['Stack','image']

        return table_cols, col_vals, band

    def _gen_cat_query(self,query_fields=None, qtype='main'):
        """
        Generate SQL Query for catalog search

        self.query is modified in place


        Args:
            query_fields (list):  Override the default list for the SQL query

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

    def get_catalog(self, query=None, query_fields=None, print_query=False, system='AB', **kwargs):
        """
        Grab a catalog of sources around the input coordinate to the search radius


        Args:
            query: Not used
            query_fields (list, optional): Over-ride list of items to query
            print_query (bool): Print the SQL query generated
            system (str): Magnitude system ['AB', 'Vega']


        Returns:
            astropy.table.Table:  Catalog of sources returned.  Includes WISE
            photometry for matched sources.
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
            main_cat = catalog_utils.ensure_empty_schema(main_cat, list(photom['VISTA'].keys()))
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


    def _select_best_img(self,imgTable,verbose,timeout=120):
        """
        Select the best band for a cutout


        Args:
            imgTable: Table of images
            verbose (bool):  Print status
            timeout (int or float):  How long to wait before timing out, in seconds


        Returns:
            HDU: header data unit for the downloaded image

        """
        row = imgTable[np.argmax(imgTable['exptime'].data.data.astype('float'))] # pick image with longest exposure time
        url = row['access_url'].decode()
        if verbose:
            print ('downloading deepest stacked image...')

        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat

    @staticmethod
    def _extract_select_map(html):
        """Extract select names and their option values from a form page."""
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
    def _extract_form_action_method(html):
        """Extract form action and method from the first form in page HTML."""
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
    def _extract_input_names(html):
        """Extract all input names from a form page."""
        input_pattern = re.compile(r'<input[^>]*name=["\']([^"\']+)["\']', re.IGNORECASE)
        return list(dict.fromkeys(input_pattern.findall(html)))

    @staticmethod
    def _extract_select_options(html, select_name):
        """Extract option values and labels for a named select element."""
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
    def _extract_fits_links(html, base_url):
        """Extract absolute FITS links from an HTML response."""
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

    def _resolve_vsa_download_link(self, session, link, timeout=120, verbose=False):
        """Resolve getImage.cgi wrapper links to direct FITS download links."""
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
    def _pick_vhs_database(database_options):
        """Pick latest VHS release value from VSA database options."""
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
    def _to_sexagesimal_strings(coord):
        """Convert ICRS coordinates to VSA-friendly sexagesimal strings."""
        ra_str = coord.ra.to_string(unit=units.hourangle, sep=':', precision=2, pad=True)
        dec_str = coord.dec.to_string(unit=units.deg, sep=':', precision=2, pad=True, alwayssign=True)
        return ra_str, dec_str

    @staticmethod
    def _choose_option(values, contains):
        """Choose first value containing token, case-insensitive."""
        token = contains.lower()
        for value in values:
            if token in value.lower():
                return value
        return None

    def _build_vsa_payload(self, html, coord, size_arcmin, band):
        """Build a permissive form payload from parsed fields and heuristics."""
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

    def _query_vsa_cutout_links(self, imsize, band, timeout=120, verbose=False):
        """Query the VSA getImage form and extract candidate FITS links."""
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
    def _select_best_vsa_link(links, band):
        """Select a deterministic best link, preferring URLs that mention band."""
        if not links:
            return None
        band_lower = band.lower()
        preferred = [link for link in links if band_lower in link.lower()]
        return preferred[0] if preferred else links[0]

    def get_image(self, imsize, band=None, timeout=120, verbose=False):
        """Retrieve a VISTA FITS image through the VSA getImage service."""
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

