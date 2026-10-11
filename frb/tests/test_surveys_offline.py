# Tests of the survey classes that do not need the network.
#  Remote queries are replaced by stubs returning small tables.

import os

import numpy as np
import pytest

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from frb.galaxies.defs import MASS_bands
from frb.surveys import catalog_utils as cu
from frb.surveys import surveycoord, dlsurvey
from frb.surveys import des, delve, nsc, vista, decals, desi, wise
from frb.surveys import sdss, hsc, euclid, skyview, twomass, panstarrs, psrcat
from frb.surveys import cluster_search

COORD = SkyCoord(10., 20., unit='deg')
RADIUS = 10 * units.arcsec


# ---------------------------------------------------------------------
# Helpers

def _dl_survey(cls, name, database, fields):
    """A DataLab survey object built without logging in to DataLab."""
    srvy = cls.__new__(cls)
    surveycoord.SurveyCoord.__init__(srvy, COORD, RADIUS)
    srvy.survey = name
    srvy.database = database
    srvy.default_query_fields = fields
    srvy.query = None
    srvy.qc_profile = 'default'
    srvy.bands = []
    srvy.svc = None
    return srvy


def _stub_dl_query(monkeypatch, table):
    """Make DL_Survey.get_catalog return (and store) a copy of ``table``."""
    def fake_get_catalog(self, query=None, query_fields=None, print_query=False,
                         timeout=120, photomdict=None):
        self.catalog = table.copy()
        return table.copy()
    monkeypatch.setattr(dlsurvey.DL_Survey, 'get_catalog', fake_get_catalog)


class _FakeTAP:
    """Stands in for a pyvo TAPService"""
    def __init__(self, table):
        self.table = table
        self.queries = []

    def run_async(self, query):
        self.queries.append(query)
        return self

    def to_table(self):
        return self.table.copy()


DL_CASES = [
    (des.DES_Survey, 'DES', 'des_dr2.main', des.photom['DES'], des.schema_dtypes['DES']),
    (delve.DELVE_Survey, 'DELVE', 'delve_dr2.objects', delve.photom['DELVE'], delve.schema_dtypes['DELVE']),
    (nsc.NSC_Survey, 'NSC', 'nsc_dr2.object', nsc.photom['NSC'], nsc.schema_dtypes['NSC']),
    (vista.VISTA_Survey, 'VISTA', 'vhs_dr5.vhs_cat_v3', vista.photom['VISTA'], vista.schema_dtypes['VISTA']),
    (decals.DECaL_Survey, 'DECaL', 'ls_dr10.tractor', decals.photom['DECaL'], decals.schema_dtypes['DECaL']),
    (desi.DESI_Survey, 'DESI', 'desi_dr1.zpix', desi.spectrom['DESI'], desi.schema_dtypes['DESI']),
]


def _check_empty(srvy, returned, schema, dtypes):
    # The stored catalog is the returned one, with the standard columns
    assert len(returned) == 0 and len(srvy.catalog) == 0
    assert srvy.catalog.colnames == returned.colnames
    for col in schema:
        assert col in returned.colnames
    for col, dtype in dtypes.items():
        assert returned[col].dtype.kind == np.dtype(dtype).kind, col
    # and it satisfies the catalog contract
    assert 'separation' in returned.colnames
    assert returned.meta['survey'] == srvy.survey
    assert returned.meta['radius'] == srvy.radius


# ---------------------------------------------------------------------
# Empty catalogs

@pytest.mark.parametrize('cls,name,database,schema,dtypes', DL_CASES)
def test_dl_empty_catalog(monkeypatch, cls, name, database, schema, dtypes):
    _stub_dl_query(monkeypatch, Table())
    srvy = _dl_survey(cls, name, database, list(schema.values()))
    _check_empty(srvy, srvy.get_catalog(), schema, dtypes)


def test_wise_empty_catalog_and_query():
    srvy = wise.WISE_Survey(COORD, RADIUS)
    srvy.service = _FakeTAP(Table())
    _check_empty(srvy, srvy.get_catalog(), wise.photom['WISE'], {})
    # The default query asks for the WISE magnitudes
    assert 'w1mag' in srvy.service.queries[-1]
    # A query given by the user is the one that is run
    srvy.get_catalog(query='SELECT ra, dec FROM allwise_p3as_psd')
    assert srvy.service.queries[-1] == 'SELECT ra, dec FROM allwise_p3as_psd'
    assert srvy.query == 'SELECT ra, dec FROM allwise_p3as_psd'


# ---------------------------------------------------------------------
# DECaLS

def _raw_decals_table():
    tab = Table()
    tab['ls_id'] = [1, 2, 3]
    tab['ra'] = [10.0001, 10.0002, 10.0003]
    tab['dec'] = [20.0001, 20.0002, 20.0003]
    tab['brickid'] = [5, 5, 5]
    tab['type'] = ['PSF', 'DEV', 'EXP']
    for band in decals.DECaL_bands:
        tab[f'mag_{band}'] = [20., 21., 22.]
        tab[f'snr_{band}'] = [50., 20., 10.]
    return tab


def test_decals_exclude_stars(monkeypatch):
    _stub_dl_query(monkeypatch, _raw_decals_table())
    fields = list(decals.photom['DECaL'].values())

    srvy = _dl_survey(decals.DECaL_Survey, 'DECaL', 'ls_dr10.tractor', fields)
    cat = srvy.get_catalog()
    assert len(cat) == 3 and 'DECaL_type' in cat.colnames

    srvy = _dl_survey(decals.DECaL_Survey, 'DECaL', 'ls_dr10.tractor', fields)
    cat = srvy.get_catalog(exclude_stars=True)
    # Point sources are removed (not kept), from the cleaned and stored catalog
    assert sorted(cat['DECaL_type']) == ['DEV', 'EXP']
    assert sorted(srvy.catalog['DECaL_type']) == ['DEV', 'EXP']
    assert 'DECaL_g' in srvy.catalog.colnames and 'separation' in srvy.catalog.colnames
    # SNR was converted to a magnitude error
    assert np.all((cat['DECaL_g_err'] > 0) & (cat['DECaL_g_err'] < 1))


# ---------------------------------------------------------------------
# Query generation

def test_vista_query_with_fields():
    srvy = _dl_survey(vista.VISTA_Survey, 'VISTA', 'vhs_dr5.vhs_cat_v3', None)
    query = srvy._gen_cat_query(query_fields=['sourceid', 'ra2000', 'dec2000'])
    assert 'sourceid' in query and 'vhs_dr5.vhs_cat_v3' in query
    # Default fields
    query = srvy._gen_cat_query()
    assert 'ypetromag' in query
    with pytest.raises(IOError):
        srvy._gen_cat_query(qtype='other')


def test_panstarrs_survey_name():
    srvy = panstarrs.Pan_STARRS_Survey(COORD, RADIUS)
    assert srvy.survey == 'Pan-STARRS'


# ---------------------------------------------------------------------
# SDSS

def test_sdss_photoz_first_row(monkeypatch):
    """The photo-z of every source is assigned, including the first one"""
    phot = Table()
    phot['ra'] = [10.0, 10.002]
    phot['dec'] = [20.0, 20.002]
    phot['objid'] = [100, 200]
    for col in ['run', 'rerun', 'camcol', 'field', 'type']:
        phot[col] = [1, 1]
    for band in sdss.SDSS_bands:
        phot[f'modelMag_{band}'] = [20., 21.]
        phot[f'modelMagErr_{band}'] = [0.1, 0.2]
        phot[f'extinction_{band}'] = [0.05, 0.05]
    photz = Table({'distance': [0.01, 0.1], 'objid': [100, 200],
                   'redshift': [0.11, 0.22], 'redshift_error': [0.01, 0.02]})

    def fake_region(coord, radius=None, timeout=None, photoobj_fields=None,
                    spectro=False, specobj_fields=None):
        return None if spectro else phot.copy()
    monkeypatch.setattr(sdss.SDSS, 'query_region', fake_region)
    monkeypatch.setattr(sdss.SDSS, 'query_sql', lambda query, timeout=None: photz.copy())

    cat = sdss.SDSS_Survey(COORD, 1 * units.arcmin).get_catalog()
    assert len(cat) == 2
    assert sorted(np.round(cat['photo_z'], 2)) == [0.11, 0.22]


# ---------------------------------------------------------------------
# HSC

def test_hsc_failed_query(monkeypatch):
    monkeypatch.setattr(hsc, 'run_query', lambda *args, **kwargs: None)
    with pytest.raises(hsc.QueryError):
        hsc.HSC_Survey(COORD, RADIUS).get_catalog()


def test_hsc_preview(monkeypatch):
    calls = []
    monkeypatch.setattr(hsc, 'getCredentials', lambda: ('user', 'password'))
    monkeypatch.setattr(hsc, '_preview', lambda credential, sql, out, release_version='pdr3': calls.append(sql))
    assert hsc.run_query('SELECT 1', preview=True) is None
    assert calls == ['SELECT 1']


# ---------------------------------------------------------------------
# Cluster catalogs

def test_cluster_empty_with_cut(monkeypatch):
    monkeypatch.setattr(cluster_search.VizierCatalogSearch, '_get_catalog',
                        lambda self, query_fields=None, **kwargs: Table(names=('ra', 'dec', 'z')))
    for cls in [cluster_search.UPClusterSZCat, cluster_search.ROSATXClusterCat,
                cluster_search.TempelClusterCat, cluster_search.RASSClusterCat,
                cluster_search.RedMapperClusterCat, cluster_search.ACTDR5ClusterCat,
                cluster_search.ERASSClusterCat]:
        srvy = cls(COORD, radius=1 * units.deg)
        cat = srvy.get_catalog(transverse_distance_cut=1 * units.Mpc)
        assert len(cat) == 0
        # The search radius keeps its units
        assert srvy.radius == 1 * units.deg
        assert cat.meta['radius'] == 1 * units.deg


def test_wen_transverse_cut(monkeypatch):
    raw = Table({'RAJ2000': [10.01, 11.0], 'DEJ2000': [20.0, 20.0],
                 'zCl': [0.1, 0.1], 'Ngal': [10, 10]})
    monkeypatch.setattr(cluster_search.VizierCatalogSearch, '_get_catalog',
                        lambda self, query_fields=None, **kwargs: raw.copy())
    srvy = cluster_search.WenGroupCat(COORD, radius=2 * units.deg)
    cat = srvy.get_catalog(transverse_distance_cut=5 * units.Mpc)
    # Only the cluster ~0.07 Mpc away in projection survives (the other is ~6 Mpc away)
    assert len(cat) == 1 and np.isclose(cat['ra'][0], 10.01)
    # Richness cut
    cat = srvy.get_catalog(richness_cut=20)
    assert len(cat) == 0


# ---------------------------------------------------------------------
# Euclid

def test_euclid_spectra_exist_no_datalinks(monkeypatch):
    monkeypatch.setattr(euclid.Euclid, 'get_datalinks', lambda ids=None: None)
    has_spec = euclid.Euclid_Survey(COORD, RADIUS).spectra_exist([1, 2])
    assert list(has_spec) == [False, False]


def test_euclid_get_spectrum(monkeypatch, tmp_path):
    srvy = euclid.Euclid_Survey(COORD, RADIUS)
    folder = str(tmp_path / 'spectra')
    seen = {}

    def fake_get_spectrum(source_id=None, output_file=None, verbose=False):
        seen['output_file'] = output_file
        return [os.path.join(folder, 'SPECTRA_RGS.fits')]
    monkeypatch.setattr(euclid.Euclid, 'get_spectrum', fake_get_spectrum)

    # A source with a spectrum: the files are returned and the folder exists
    monkeypatch.setattr(srvy, 'spectra_exist', lambda ids: np.array([True]))
    files = srvy.get_spectrum(123, output_folder=folder)
    assert files == [os.path.join(folder, 'SPECTRA_RGS.fits')]
    assert os.path.isdir(folder) and seen['output_file'].startswith(folder)

    # No spectrum, or a failure: an empty list
    monkeypatch.setattr(srvy, 'spectra_exist', lambda ids: np.array([False]))
    assert srvy.get_spectrum(123, output_folder=folder) == []

    def broken(ids):
        raise RuntimeError('archive down')
    monkeypatch.setattr(srvy, 'spectra_exist', broken)
    assert srvy.get_spectrum(123, output_folder=folder) == []


# ---------------------------------------------------------------------
# NSC images

class _FakeSIA:
    """Stands in for a pyvo SIAService"""
    def __init__(self, table):
        self.table = table
        self.sizes = []

    def search(self, coord, size, verbosity=2):
        self.sizes.append(size)
        return self

    def to_table(self):
        return self.table.copy()


def _fake_cutout(nx, ny, pixscale=0.27):
    """A cutout of nx x ny pixels with a simple WCS"""
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [COORD.ra.deg, COORD.dec.deg]
    w.wcs.crpix = [nx / 2, ny / 2]
    w.wcs.cdelt = [-pixscale / 3600, pixscale / 3600]
    return fits.PrimaryHDU(data=np.zeros((ny, nx)), header=w.to_header())


def _nsc_images():
    # CCD images: centre offsets (arcmin) from the target, depth, band, type
    rows = [  # name,   offset, magzero, band, proctype,    prodtype
        ('deep_far',     6.0,  '30.0',  'g',  'Resampled', 'image'),
        ('shallow_near', 1.0,  '28.0',  'g',  'Resampled', 'image'),
        ('deep_near',    2.0,  '29.0',  'g',  'Resampled', 'image'),
        ('other_band',   1.0,  '31.0',  'r',  'Resampled', 'image'),
        ('instcal',      1.0,  '31.0',  'g',  'InstCal',   'image'),
        ('weight',       1.0,  '31.0',  'g',  'Resampled', 'wtmap'),
    ]
    tab = Table()
    tab['s_ra'] = [COORD.ra.deg] * len(rows)
    tab['s_dec'] = [COORD.dec.deg + row[1] / 60 for row in rows]
    tab['magzero'] = np.ma.MaskedArray([row[2] for row in rows], mask=False)
    tab['obs_bandpass'] = [row[3] for row in rows]
    tab['proctype'] = [row[4] for row in rows]
    tab['prodtype'] = [row[5] for row in rows]
    tab['access_url'] = [f'https://example.org/cutout?siaRef={row[0]}&POS=10.0,20.0&SIZE=0.2,0.2'
                         for row in rows]
    return tab


def _nsc_survey(monkeypatch, cutouts):
    """NSC survey whose image service and downloads are stubbed.

    ``cutouts`` maps the name of an image to the (nx, ny) of the cutout it returns.
    """
    srvy = _dl_survey(nsc.NSC_Survey, 'NSC', 'nsc_dr2.object', None)
    srvy.bands = nsc.NSC_bands
    srvy.svc = _FakeSIA(_nsc_images())
    srvy.urls = []

    def fake_download(url, timeout=120):
        srvy.urls.append(url)
        name = url.split('siaRef=')[1].split('&')[0]
        return _fake_cutout(*cutouts[name])
    monkeypatch.setattr(srvy, '_download_cutout', fake_download)
    return srvy


def test_nsc_image_any_size(monkeypatch):
    """The cutout is requested at the size asked for, whatever it is"""
    for imsize in [5 * units.arcsec, 10 * units.arcsec, 2 * units.arcmin]:
        npix = int(round(imsize.to('arcsec').value / 0.27))
        srvy = _nsc_survey(monkeypatch, {'deep_near': (npix, npix), 'shallow_near': (npix, npix),
                                         'deep_far': (npix, npix)})
        hdu = srvy.get_image(imsize, 'g')
        assert isinstance(hdu, fits.PrimaryHDU) and hdu.data.shape == (npix, npix)
        # The search region is wide even for a small image...
        assert srvy.svc.sizes[-1] >= 0.2 * units.deg
        # ...but the cutout has the requested size
        size = imsize.to('deg').value
        assert f'SIZE={size:.6f},{size:.6f}' in srvy.urls[-1]
        # The deepest image sure to contain the cutout is used, in one download:
        # not the deeper one whose centre is too far, nor other bands or products
        assert len(srvy.urls) == 1 and 'siaRef=deep_near' in srvy.urls[0]


def test_nsc_image_truncated_cutout(monkeypatch):
    """A cutout truncated by the edge of the CCD is replaced by a complete one"""
    srvy = _nsc_survey(monkeypatch, {'deep_near': (12, 37), 'shallow_near': (37, 37),
                                     'deep_far': (37, 37)})
    hdu = srvy.get_image(10 * units.arcsec, 'g')
    assert hdu.data.shape == (37, 37)
    assert [url.split('siaRef=')[1].split('&')[0] for url in srvy.urls] == ['deep_near', 'shallow_near']

    # Nothing complete: the best coverage is returned with a warning
    srvy = _nsc_survey(monkeypatch, {'deep_near': (12, 37), 'shallow_near': (37, 20),
                                     'deep_far': (5, 5)})
    with pytest.warns(RuntimeWarning, match='covers only'):
        hdu = srvy.get_image(10 * units.arcsec, 'g')
    assert hdu.data.shape == (20, 37)


def test_nsc_image_band_handling(monkeypatch):
    srvy = _nsc_survey(monkeypatch, {'other_band': (37, 37)})
    # Default band is r; the band is case-insensitive
    assert srvy.get_image(10 * units.arcsec).data.shape == (37, 37)
    assert 'siaRef=other_band' in srvy.urls[-1]
    assert srvy.get_image(10 * units.arcsec, 'R') is not None
    # No image in this band
    assert srvy.get_image(10 * units.arcsec, 'z') is None
    with pytest.raises(TypeError):
        srvy.get_image(10 * units.arcsec, 'q')
    # The deprecated alias goes through get_image
    with pytest.warns(DeprecationWarning):
        data, hdr = srvy.get_cutout(10 * units.arcsec, 'r')
    assert data.shape == (37, 37) and isinstance(hdr, fits.Header)

# ---------------------------------------------------------------------
# SkyView

def test_skyview_gleam_band(monkeypatch):
    srvy = skyview.SkyView_Survey(COORD, RADIUS, 'gleam')
    names = []
    monkeypatch.setattr(srvy, '_skyview_fetch',
                        lambda name, imsize, pixels=None: names.append(name))
    srvy.get_image(imsize=1 * units.deg, band='72-103 MHz')
    srvy.get_image(imsize=1 * units.deg)
    assert names == ['GLEAM 72-103 MHz', 'GLEAM 170-231 MHz']


def test_skyview_at_least_one_pixel(monkeypatch):
    srvy = skyview.SkyView_Survey(COORD, RADIUS, 'nvss')
    seen = {}

    def fake_get_images(position=None, survey=None, radius=None, pixels=None):
        seen['pixels'] = pixels
        return []
    monkeypatch.setattr(skyview.SkyView, 'get_images', fake_get_images)
    # 1 arcsec is far below the 15 arcsec/pixel scale of NVSS
    with pytest.warns(UserWarning):
        assert srvy._skyview_fetch('NVSS', 1 * units.arcsec) is None
    assert seen['pixels'] == '1'


# ---------------------------------------------------------------------
# 2MASS, Pan-STARRS, PSRCAT

def test_twomass_data_model_untouched(monkeypatch):
    psc = Table({'designation': ['00400000+2000000'], 'ra': [10.0001], 'dec': [20.0001]})
    for band in MASS_bands:
        psc[f'{band}_m'] = [15.]
        psc[f'{band}_msigcom'] = [0.05]

    def fake_region(coord, radius=None, spatial=None, catalog=None):
        return Table() if catalog == 'fp_xsc' else psc.copy()
    monkeypatch.setattr(twomass.Irsa, 'query_region', fake_region)

    before = dict(twomass.photom['2MASS'])
    cat = twomass.TwoMASS_Survey(COORD, RADIUS).get_catalog()
    # The point source catalog was used and its columns were renamed
    assert len(cat) == 1 and '2MASS_j' in cat.colnames and '2MASS_j_err' in cat.colnames
    # without changing the module-level data model
    assert twomass.photom['2MASS'] == before


def test_ps1_metadata_cache_not_writable(monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return [{'name': 'objID', 'datatype': 'bigint', 'description': 'ID'}]

    def no_write(self, *args, **kwargs):
        raise PermissionError('read-only installation')
    monkeypatch.setattr(panstarrs.requests, 'get', lambda url: FakeResponse())
    monkeypatch.setattr(panstarrs.Table, 'write', no_write)
    tab = panstarrs._ps1metadata()
    assert list(tab['name']) == ['objID']


def test_psrcat_missing_dependency(monkeypatch):
    monkeypatch.setattr(psrcat, 'pio', None)
    with pytest.raises(ImportError):
        psrcat.PSRCAT_Survey(COORD, RADIUS).get_catalog()


# ---------------------------------------------------------------------
# catalog_utils

def test_convert_mags_missing_error_column():
    """Each magnitude is paired with its own error, even if another has none"""
    tab = Table({'DES_g': [20.], 'DES_r': [21.], 'DES_r_err': [0.1]})
    flux = cu.convert_mags_to_flux(tab, 'mJy')
    assert np.isclose(flux['DES_g'][0], 3630.7805e3 * 10**(-20 / 2.5))
    assert np.isclose(flux['DES_r'][0], 3630.7805e3 * 10**(-21 / 2.5))
    assert np.isclose(flux['DES_r_err'][0], np.log(10) / 2.5 * flux['DES_r'][0] * 0.1)


def test_zero_errors_are_bad():
    # Flux conversion: the flux is kept, its error is flagged
    tab = Table({'DES_r': [20., 21.], 'DES_r_err': [0., 0.1]})
    flux = cu.convert_mags_to_flux(tab, 'mJy')
    assert flux['DES_r'][0] > 0 and flux['DES_r_err'][0] == -99.
    assert flux['DES_r_err'][1] > 0
    # Survey cleaning: the same
    tab = Table({'ra': [1., 2.], 'dec': [1., 2.], 'mag': [20., 21.], 'err': [0., 0.1]})
    out = cu.clean_cat(tab, {'T_g': 'mag', 'T_g_err': 'err', 'ra': 'ra', 'dec': 'dec'},
                       mask_photometry=True)
    assert out['T_g'][0] == 20. and out['T_g_err'][0] == -99.
    assert out['T_g_err'][1] == 0.1


def test_summarize_catalog_closest():
    cat = Table({'ra': [10.0001, 10.001], 'dec': [20., 20.], 'flux': [1., 10.]})
    cat.meta['survey'] = 'TEST'
    summary = cu.summarize_catalog({'coord': COORD}, cat, 10 * units.arcsec, 'flux', False)
    assert 'brightest source has flux of 10.00' in summary[1]
    # The closest source is the faint one
    assert 'flux of 1.00' in summary[2]


def test_xmatch_catalogs_checks_both():
    tab = Table({'ra': [1.], 'dec': [1.]})
    with pytest.raises(AssertionError):
        cu.xmatch_catalogs(tab, 'not a table')
