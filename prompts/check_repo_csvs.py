""" Consistency checks for FRBs_base.csv, public_hosts.csv and the
literature tables (add_frbs.md, prompt 12).  Read-only.

Run from the repo root:  python prompts/check_repo_csvs.py
"""
import glob
import json
import os
import re

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

FRB_PATH = 'frb/data/FRBs/'
GAL_PATH = 'frb/data/Galaxies/'
LIT_PATH = GAL_PATH + 'Literature/'

frbs = pd.read_csv(FRB_PATH + 'FRBs_base.csv')
hosts = pd.read_csv(GAL_PATH + 'public_hosts.csv', dtype={'FRB': str})

# Basic content
print('Shapes:', frbs.shape, hosts.shape)
print('Duplicate names:', frbs.Name[frbs.Name.duplicated()].tolist(),
      hosts.FRB[hosts.FRB.duplicated()].tolist())
print('Non-standard columns:', [c for c in frbs.columns if c.endswith('_new') or c.startswith('Diff')])
for col in ['ra', 'dec', 'ee_a', 'ee_b']:
    print(f'FRBs missing {col}:', frbs.loc[frbs[col].isna(), 'Name'].tolist())
print('Hosts missing Coord:', hosts.loc[hosts.Coord.isna(), 'FRB'].tolist())
print('Hosts missing References:', hosts.loc[hosts.References.isna(), 'FRB'].tolist())

# FRB JSONs vs FRBs_base.csv
for _, row in frbs.iterrows():
    jfile = FRB_PATH + f'{row.Name}.json'
    if not os.path.isfile(jfile):
        print('No FRB JSON:', row.Name)
        continue
    d = json.load(open(jfile))
    zb = row.z if np.isfinite(row.z) else None
    z = d.get('z')
    if (z is None) != (zb is None) or (z is not None and abs(z - zb) > 1e-6):
        print('FRB JSON z differs:', row.Name, zb, z)
    if abs(d['ra'] - row.ra) > 1e-5 or abs(d['dec'] - row.dec) > 1e-5:
        print('FRB JSON position differs:', row.Name)
json_names = [os.path.basename(f)[:-5] for f in glob.glob(FRB_PATH + 'FRB*.json')]
print('FRB JSONs not in FRBs_base.csv:', sorted(set(json_names) - set(frbs.Name)))

# Host JSONs vs public_hosts.csv
host_dirs = [os.path.basename(os.path.dirname(f))
             for f in glob.glob(GAL_PATH + '*/FRB*_host.json')]
print('Host dirs not in public_hosts.csv:', sorted(set(host_dirs) - set(hosts.FRB)))
print('Hosts not in FRBs_base.csv:', [n for n in hosts.FRB if 'FRB' + n not in set(frbs.Name)])
fidx = frbs.set_index('Name')
for _, row in hosts.iterrows():
    if not isinstance(row.Coord, str):
        continue
    if not re.match(r'^\d{2}h\d{2}m[\d.]+s [+-]\d{2}d\d{2}m[\d.]+s$', row.Coord):
        print('Non-standard Coord:', row.FRB, repr(row.Coord))
    coord = SkyCoord(row.Coord, frame='icrs')
    jfile = GAL_PATH + f'{row.FRB}/FRB{row.FRB}_host.json'
    if not os.path.isfile(jfile):
        print('No host JSON:', row.FRB)
    else:
        d = json.load(open(jfile))
        sep = coord.separation(SkyCoord(d['ra'], d['dec'], unit='deg')).arcsec
        if sep > 0.25:
            print(f'Host JSON position off by {sep:.2f}":', row.FRB)
    # Host vs FRB position (catches a Coord that parses to the wrong place)
    name = 'FRB' + row.FRB
    if name in fidx.index and np.isfinite(fidx.loc[name, 'ra']):
        sep = coord.separation(SkyCoord(fidx.loc[name, 'ra'], fidx.loc[name, 'dec'],
                                        unit='deg')).arcsec
        if sep > 30 and sep > 5 * fidx.loc[name, 'ee_a']:
            print(f'Host is {sep:.0f}" from the FRB (ee_a = {fidx.loc[name, "ee_a"]}"):', row.FRB)
    # z
    if name in fidx.index:
        zf, zh = fidx.loc[name, 'z'], row.z
        if not np.isclose(np.nan_to_num(zf, nan=-1), np.nan_to_num(zh, nan=-1), atol=1e-5):
            print('z differs (FRBs_base, public_hosts):', row.FRB, zf, zh)
        pf = fidx.loc[name, 'P(O|x)']
        if np.isfinite(pf) and not np.isclose(pf, row.P_Ox, atol=1e-3):
            print('P(O|x) differs:', row.FRB, pf, row.P_Ox)

# Literature tables
good_hosts = hosts.dropna(subset=['Coord']).reset_index()
hcoord = SkyCoord(list(good_hosts.Coord), frame='icrs')
refs = pd.read_csv(LIT_PATH + 'all_refs.csv', comment='#')
for _, entry in refs.iterrows():
    if not os.path.isfile(LIT_PATH + entry.Table):
        print('Missing literature table:', entry.Table)
        continue
    if entry.Format == 'csv':
        tbl = Table.from_pandas(pd.read_csv(LIT_PATH + entry.Table))
    else:
        tbl = Table.read(LIT_PATH + entry.Table, format=entry.Format)
    idx, sep, _ = SkyCoord(tbl['ra'], tbl['dec'], unit='deg').match_to_catalog_sky(hcoord)
    for ii in np.where(sep > 1 * u.arcsec)[0]:
        print(f'{entry.Table}: row {ii} has no host within 1" (nearest {good_hosts.FRB[idx[ii]]})')
    for col in tbl.colnames:
        if col.endswith('_loerr'):
            vals = np.array(tbl[col], dtype=float)
            nneg = np.sum(np.isfinite(vals) & (vals < 0) & ~np.isin(vals, [-999, -998]))
            if nneg:
                print(f'{entry.Table}: {nneg} negative {col}')
