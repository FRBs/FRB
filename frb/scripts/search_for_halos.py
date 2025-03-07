"""
Script to search any given sightline for all the halo intersections within
the NEDLVS and Tully 2015 catalogs.
"""
from frb.halos.models import ModifiedNFW, ICM, halomass_from_stellarmass_kravtsov #, YF17
from frb.defs import frb_cosmo
from frb.surveys.catalog_utils import xmatch_catalogs

from astropy import units as u
from astropy.table import Table, setdiff
from astropy.coordinates import SkyCoord
from astropy.cosmology import z_at_value

import numpy as np

def parser(options=None):
  import argparse
  # Parse
  parser = argparse.ArgumentParser(description='Script to search for halos along a given sightline')
  parser.add_argument("frb_name", type=str, help="Name of the FRB")
  parser.add_argument("frb_coord", type=str, help="Coordinates of the FRB. e.g. 122.223,-23.2322")
  parser.add_argument("frb_z", type=float, help="Redshift of the FRB")
  parser.add_argument("--rmax", type=float, default=1., help="Radial extent of halo model in units of r_vir [default=1.].")

  if options is None:
      pargs = parser.parse_args()
  else:
      pargs = parser.parse_args(options)
  return pargs


def halo_dm(z, offset, log_mhalo, rmax=1, fhot=1, is_grp=False):
  if log_mhalo<14: # Low mass halos
      if offset>1300*u.kpc: # Don't bother instantiating the model if the distance is too far
          return 0, -99*u.kpc
      mnfw = ModifiedNFW(log_Mhalo = log_mhalo, alpha = 2, y0 = 2, z = z, f_hot=fhot)
      #mnfw = YF17(log_Mhalo = log_mhalo, z = z)
  elif log_mhalo>=14 & is_grp: # Clusters
      mnfw = ICM(log_mhalo, f_hot=fhot, z=z)
      import pdb; pdb.set_trace()
#         if distance<mnfw.r200: # No local group
#             return 0

  dm_halo, rvir = mnfw.Ne_Rperp(offset, rmax=rmax).to('pc/cm**3').value/(1+z), mnfw.r200
  return dm_halo, rvir
    
def lvs_avg_dm_halos(frb_name, frb_coord, frb_z, nedlvs_tab, tully_clusters, rmax=1):
  """
  Estimate the lcontribution of the local volume halos available
  within the NEDLVS catalog. Produces estimates where masses are available
  """
  print(f"Working on {frb_name}")
  
  # Conditions for closeness
  #import pdb; pdb.set_trace()
  valid_distances = nedlvs_tab['DistMpc']>0 # Weed out weird ones
  distance_cut = nedlvs_tab['DistMpc']<frb_cosmo.luminosity_distance(frb_z).to('Mpc').value #Only need foreground objects
  #local_grp_cut = nedlvs_tab['DistMpc']>2 # Exclude local group
  phys_sep_cut = nedlvs_tab['phys_sep']<1*u.Mpc # Impact param within 1 Mpc
  ang_sep_cut = nedlvs_tab['ang_sep']<90*u.deg # Make sure the earth is not between the FRB and the galaxy
  is_nearby_fg = valid_distances&distance_cut & phys_sep_cut & ang_sep_cut
  
  # Get culled table
  close_by = nedlvs_tab[is_nearby_fg]
  print(f"{len(close_by)} objects are identified in the foreground of {frb_name}")
  
  # Have mass estimates?
  mass_cut = (close_by['Mstar']>0) & (close_by['Mstar']<1e12) # Weird to have galaxies above this value
  close_by_withmass = close_by[mass_cut]
  print(f"{len(close_by_withmass)} objects already have masses")
      
  # If no masses, then inform the user.
  if len(close_by_withmass)==0:
      print("No local volume halos with mass found in the FRB f/g.")
  
  # Get DM estimates
  logmstar = np.log10(close_by_withmass["Mstar"]/1.7) # Start with Mstar. Reduce by 0.3 dex for Salpeter to Chabrier IMF

  logmstar[logmstar>10.66] += 0.2 # Add 0.2 dex to "correct" for missing light if using the Kravtsov SHMR. THIS IS AD HOC.
  
  logmhalo = halomass_from_stellarmass_kravtsov(logmstar) # Convert to Mhalo
  close_by_withmass['log_mhalo'] = logmhalo 
  close_by_withmass = close_by_withmass.filled(-99.)
  
  rvirs = []
  halo_dms = []
  # Loop through the objects and populate the DMs in the table.
  for num, (z, lmh, offset, distance) in enumerate(zip(close_by_withmass['z'], logmhalo, close_by_withmass['phys_sep'], close_by_withmass['DistMpc'])):
      if z<0 and close_by_withmass[num]['DistMpc_method'] == 'zIndependent':
          z = 0
      dmh, rvir = halo_dm(z=z,
                          offset=offset*u.Mpc,
                          log_mhalo=lmh, rmax=rmax)
      if np.abs(frb_z-z)<=2e-3:
          #dmh /= 2. # In case they're at the same z, then just assume only half the halo is intersected
          dmh = 0.0 # Strictly fg
      halo_dms.append(dmh)
      rvirs.append(rvir.to('kpc'))
  close_by_withmass['DM_halo'] = halo_dms
  close_by_withmass['rvir'] = rvirs
  close_by_withmass.sort("DM_halo")
  
  # Are there any group members in the field?
  # Compute transverse distances and initialize columns
  match_grps = tully_clusters[np.isin(tully_clusters['PGC'], np.unique(tully_clusters['PGC1']))] # This gets the BCGs only.
  match_grps['coord'] = SkyCoord(ra=match_grps['ra'], dec=match_grps['dec'], unit='deg')
  match_grps['ang_sep'] = frb_coord.separation(match_grps['coord']).to('arcmin')
  match_grps['phys_sep'] = match_grps['Dist']*np.sin(match_grps['ang_sep'].to('rad').value)
  match_grps = match_grps[match_grps['phys_sep']<5*u.Mpc] # Only consider groups within 5 Mpc
  
  # Cross match with the tully cluster catalog
  if len(match_grps)>0:
      use_grps = True
      match_grps['DM_halo'] = 0.0
      match_grps['rvir'] = 0.0*u.kpc
      match_lvs, _ = xmatch_catalogs(close_by_withmass, tully_clusters, skydist=1*u.arcsec)

      if len(match_lvs)>0:
          # remove matched objects from the table
          close_by_withmass = setdiff(close_by_withmass, match_lvs,keys="objname")
  else:
      use_grps = False

  print(len(match_grps), " groups found within 5Mpc.")
  # Loop through groups
  if use_grps:
      for central_entry in match_grps:
          log_grp_mass = np.log10(central_entry['Mlum'])+12. # Mlum is in Tera Msun units.
          z_grp = z_at_value(frb_cosmo.luminosity_distance, central_entry['Dist']*u.Mpc)
          dm_halo, rvir = halo_dm(z=z_grp,
                                  offset=central_entry['phys_sep']*u.Mpc,
                                  log_mhalo=log_grp_mass, is_grp=True, rmax=rmax)
          if np.abs(z_grp-frb_z)<=2e-3:
              #dm_halo /= 2.
              dm_halo=0.0 # Strictly foreground sources
          central_entry['DM_halo'] = dm_halo
          central_entry['rvir'] = rvir.to('kpc').value
  else:
      match_grps = Table()

  
  # return average dm_halos
  mean_dm_halos_lvs = np.sum(close_by_withmass["DM_halo"])
  if use_grps:
      mean_grp_dm = np.sum(match_grps['DM_halo'])
  else:
      mean_grp_dm = 0.0
  photo_z_halos = (close_by_withmass['z_tech']=='PHOT')&(close_by_withmass['DistMpc_method']=='Redshift')
  mean_dm_halos_lvs_phot = np.sum(close_by_withmass['DM_halo'][photo_z_halos])
  print("<DM_halos> = ", mean_dm_halos_lvs)
  print("<DM_halos>_photoz = ", mean_dm_halos_lvs_phot)
  print("<DM_halos>_grp = ", mean_grp_dm)
  print("-"*20)
  return mean_dm_halos_lvs, mean_dm_halos_lvs_phot, mean_grp_dm, close_by_withmass, match_grps

def main(pargs):

  from frb.surveys.cluster_search import TullyGroupCat
  from frb.surveys.nedlvs import NEDLVS
  
  # First read in the coords
  ra, dec = pargs.frb_coord.split(',')
  ra = float(ra)
  dec = float(dec)
  frb_coord = SkyCoord(ra=ra, dec=dec, unit="deg")
  frb_name = pargs.frb_name
  frb_z = pargs.frb_z
  rmax = pargs.rmax

  nedlvs_tab = NEDLVS(frb_coord, radius=90*u.deg,).get_catalog(z_lim=frb_z, impact_par_lim=1*u.Mpc)

  tully_clusters = TullyGroupCat(frb_coord, radius=90*u.deg).get_catalog()

  mean_dm_halos_lvs, mean_dm_halos_lvs_phot, mean_grp_dm, close_by_withmass, match_grps = lvs_avg_dm_halos(frb_name, 
                                                                                                           frb_coord, 
                                                                                                           frb_z, 
                                                                                                           nedlvs_tab, 
                                                                                                           tully_clusters, 
                                                                                                           rmax=rmax)

  close_by_withmass.write(f"{frb_name}_nedvls_with_mass.fits", overwrite=True)
  if len(match_grps)>0:
    match_grps.write(f"{frb_name}_tully_grps_ICM_fgonly.fits", overwrite=True)
  