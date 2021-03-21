""" Module for methods related to building FRB Assocation sandboxes
Warning:  These methods are effectively deprecated. """

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units

from frb import frb

from IPython import embed

def theta_rcore(theta, rstate, ntheta):
    rand = rstate.uniform(size=ntheta)
    dtheta = np.exp(np.log(theta['core']) + rand) - theta['core']
    embed(header='16 of sandbox')
    pa = rstate.uniform(size=ntheta, low=0., high=360.)
    # Return
    return dtheta*units.arcsec, pa*units.deg

def theta_uniform(theta, rstate, ntheta):
    dtheta = rstate.uniform(size=ntheta, low=0., high=theta['max'])
    pa = rstate.uniform(size=ntheta, low=0., high=360.)
    # Return
    return dtheta*units.arcsec, pa*units.deg

def offset_frb(sigR, rstate, nFRB, max_sig=5.):
    roff = rstate.normal(scale=sigR, size=nFRB)
    pa = rstate.uniform(size=nFRB, low=0., high=360.)
    # Turn into ra, dec
    roff = np.abs(np.minimum(roff, max_sig*sigR))  # Note this includes negative values
    # Return
    return roff*units.arcsec, pa*units.deg


def build(source_tbl, field, phot_col, zmnx, Lmnx, theta, sigR, nSand=100,
          seed=12345, edge_buff=20*units.arcsec, rmag_mnx=(19., 22.)):

    # Random numbers
    rstate = np.random.RandomState(seed)

    # Choose random redshifts
    zval = rstate.uniform(size=nSand, low=zmnx[0], high=zmnx[1])

    # Convert z to r_mag from Lmnx.  Cheating for now
    rmags = rstate.uniform(size=nSand, low=rmag_mnx[0], high=rmag_mnx[1])

    # Sandbox edge
    field_coord = SkyCoord(ra=field[0], dec=field[1], unit='deg')
    obj_coord = SkyCoord(ra=source_tbl['ra'], dec=source_tbl['dec'], unit='deg')
    obj_sep = field_coord.separation(obj_coord)
    ok_sep = obj_sep < (field[2]*units.deg - edge_buff)

    # Collect the objects
    chosen = np.zeros(len(source_tbl), dtype=bool)
    delta_r = 0.2  # mag

    objs = []
    for kk in range(nSand):
        rmag = rmags[kk]
        # Pick one
        options = (np.abs(source_tbl[phot_col] - rmag) < delta_r) & np.logical_not(chosen) & ok_sep
        noptions = np.sum(options)
        if noptions < 10:
            embed(header='29 of sandbox')

        # Random select
        rand = int(rstate.randint(low=0, high=noptions-1, size=1))
        iobj = np.where(options)[0][rand]
        objs.append(iobj)
        chosen[iobj] = True
    objs = np.array(objs)

    # Offset
    if theta['method'] == 'rcore':
        thetas_off, thetas_pa = theta_rcore(theta, rstate, nSand)
    elif theta['method'] == 'uniform':
        thetas_off, thetas_pa = theta_uniform(theta, rstate, nSand)
    else:
        embed(header='45')

    # FRB offset
    frb_off, frb_pa = offset_frb(sigR, rstate, nSand)

    # Apply
    frb_coords = obj_coord[objs]
    obs_coords = []
    for kk, frb_coord in enumerate(frb_coords):
        obs_coord = frb_coord.directional_offset_by(thetas_pa[kk], thetas_off[kk])
        obs_coord = obs_coord.directional_offset_by(frb_pa[kk], frb_off[kk])
        obs_coords.append(obs_coord)
    frb_coords = SkyCoord(obs_coords)

    # Final table
    frb_tbl = Table()
    frb_tbl[phot_col] = rmags
    frb_tbl['ra'] = frb_coords.ra.value
    frb_tbl['dec'] = frb_coords.dec.value
    frb_tbl['iobj'] = objs
    frb_tbl['obj_ra'] = obj_coord[objs].ra.value
    frb_tbl['obj_dec'] = obj_coord[objs].dec.value
    frb_tbl['theta'] = thetas_off.value

    # Return
    return frb_tbl

# Command line execution
if __name__ == '__main__':
    if False:  # Bright
        # Test
        source_tbl = Table.read('dev/tst_DES_180924.fits')
        field = source_tbl.meta['RA'], source_tbl.meta['DEC'], source_tbl.meta['RSEARCH']
        #theta = dict(method='rcore', max=2., core=0.1)
        theta = dict(method='uniform', max=2.)
        frb_tbl = build(source_tbl, field, 'DES_r', (0.2, 0.4), (0.1, 1.), theta, 0.25,
                        rmag_mnx = (19., 22.))
        frb_tbl.write('dev/tst_FRB_180924.fits', overwrite=True)

    if False:  # Faint
        # Test
        source_tbl = Table.read('dev/tst_DES_180924.fits')
        field = source_tbl.meta['RA'], source_tbl.meta['DEC'], source_tbl.meta['RSEARCH']
        theta = dict(method='uniform', max=2.)
        frb_tbl = build(source_tbl, field, 'DES_r', (0.2, 0.4), (0.1, 1.), theta, 0.25,
                        rmag_mnx=(21., 24.))
        frb_tbl.write('dev/tst_FRB_180924_faint.fits', overwrite=True)

    if False:  # Faint + theta_max = 3.5
        # Test
        source_tbl = Table.read('dev/tst_DES_180924.fits')
        field = source_tbl.meta['RA'], source_tbl.meta['DEC'], source_tbl.meta['RSEARCH']
        theta = dict(method='uniform', max=3.5)
        frb_tbl = build(source_tbl, field, 'DES_r', (0.2, 0.4), (0.1, 1.), theta, 0.25,
                        rmag_mnx=(21., 24.))
        frb_tbl.write('dev/tst_FRB_180924_faint_theta3.5.fits', overwrite=True)

    if True:  # Faint + sigR=0.75"
        # Test
        source_tbl = Table.read('dev/tst_DES_180924.fits')
        field = source_tbl.meta['RA'], source_tbl.meta['DEC'], source_tbl.meta['RSEARCH']
        theta = dict(method='uniform', max=2.)
        frb_tbl = build(source_tbl, field, 'DES_r', (0.2, 0.4), (0.1, 1.), theta, 0.75,
                        rmag_mnx=(21., 24.))
        frb_tbl.write('dev/tst_FRB_180924_faint_sigR0.75.fits', overwrite=True)
