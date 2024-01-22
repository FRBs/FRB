import copy
import csv
import math
import os

import fitsio
import matplotlib.pyplot as plt
import numpy as np
import sep
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from scipy import signal


def get_star_norm(size, FWHM):
	"""This function gets a normalized gaussian matrix"""

	std = FWHM/2.35482
	get_star_1d = signal.gaussian(size, std=std).reshape(size, 1)
	get_star_2d = np.outer(get_star_1d, get_star_1d)
	return get_star_2d

def get_star(star_norm, flux):
	"""This function gets a gaussian matrix with a determined flux"""

	# valor A que debo multiplicar a la matriz para que la suma de elementos sea igual al flujo:
	A = flux / np.sum(star_norm)
	star = star_norm * A
	return star

def get_mag_rec(image_diff, CC, circle_radio, x, y):
	"""This function gets the magnitude of a star at position x, y of a differentiated image"""

	data = fitsio.read(image_diff)
    # sustraccion background
	bkg = sep.Background(data)     # estima el background
	data_sub = data - bkg    # resta el background estimado
	objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)    # detecta objetos
	flux_s, fluxerr_s, flag_s = sep.sum_circle(data_sub, x, y, circle_radio, err=bkg.globalrms, gain=1.0)

	flux_s = np.where(flux_s>0, flux_s, 1)
	mas = flux_s + fluxerr_s
	menos = flux_s - fluxerr_s
	mas = np.where(mas>0, mas, 1)
	menos = np.where(menos>0, menos,1)
	mag = CC -2.5 * math.log10(flux_s)
	mag_m = CC -2.5 * math.log10(mas)    # porque mas flujo es menos magnitude
	mag_p = CC -2.5 * math.log10(menos)

	err_mag_m = mag - mag_m
	err_mag_p = mag_p - mag

	SN = flux_s/fluxerr_s
	
	return (mag, err_mag_p, err_mag_m, SN)

def mag_rec(image, skymapper, out_image_star, filtro, mag, x, y, h_size, FWHM, circle_radio, image_template, image_diff):
	"""This function makes an artificial star of a given magnitude at the x, y position of an image, differentiates this 
	image (ToO) with a template from a different epoch and calculates the recovered magnitude of the artificial star in 
	the difference image ."""

	CC = get_CC(image, skymapper, filtro)
	flux = get_flux(mag, CC)
	star_norm = get_star_norm(2*h_size, FWHM)
	star = get_star(star_norm, flux)

	zeros = np.zeros((4096, 4096))
	y_menos = y - h_size
	y_mas = y + h_size
	x_menos = x - h_size
	x_mas = x + h_size
	zeros[y_menos:y_mas, x_menos:x_mas] = star

	data = fitsio.read(image)
	copy_data = copy.deepcopy(data)
	copy_data_star = copy_data + zeros
	hdu = fits.PrimaryHDU(copy_data_star)
	hdulist = fits.HDUList([hdu])
	hdulist.writeto(out_image_star)

	image_ToO = out_image_star

	os.system("hotpants "\
	 			"-inim {} "\
 				"-tmplim {} "\
 				"-outim {} "\
 				"-c i "\
 				"-tu 50000 "\
 				"-iu 50000 "\
 				"-tl -500 "\
 				"-il -500 "\
 				"-tr 9.5 "\
 				"-ir 9.5".format(image_ToO, image_template, image_diff))

	mag_tup = get_mag_rec(image_diff, CC, circle_radio, x, y) 
	return mag_tup










def create_median(cube):
	cube_array = fitsio.read(cube)
	median = np.median(cube_array, axis=0)
	return median

def get_fits_images(array_image, suffix):
	hdu = fits.PrimaryHDU(array_image)
	hdulist = fits.HDUList([hdu])                # hdul = header data unit list
	fits_path = '/home/consuelo/projects/FRBs/alopeke/images/image' + '_' + suffix + '.fits'
	hdulist.writeto(fits_path, overwrite = True)
	return fits_path

def get_diff(image, template, image_diff, convolve_with='i'):
	os.system("hotpants "\
 				"-inim {} "\
				"-tmplim {} "\
				"-outim {} "\
				"-c {} "\
				"-tu 60000 "\
				"-iu 60000 "\
				"-tl -10 "\
				"-il -10 ".format(image, template, image_diff, convolve_with))


def search_obj(image, plot = False):
	img = fitsio.read(image)
	bkg = sep.Background(img)
	img_sub = img - bkg
	
	objects = sep.extract(img_sub, 3, err=bkg.globalrms)

	if plot:
		fig, ax = plt.subplots()
		m = np.mean(img_sub)
		s = np.std(img_sub)
		im = ax.imshow(img_sub, interpolation='nearest', cmap='gray', vmin = m-s, vmax = m+s, origin='lower')

		# plot an ellipse for each object
		for i in range(len(objects)):
		    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
		                width=6*objects['a'][i],
		                height=6*objects['b'][i],
		                angle=objects['theta'][i] * 180. / np.pi)
		    e.set_facecolor('none')
		    e.set_edgecolor('red')
		    ax.add_artist(e)

		plt.show()

	return objects


def check_objects_in_region(image, objects, xmin, xmax, ymin, ymax):
	all_objects = objects
	bool_aux = False
	for x,y in zip(all_objects['x'],all_objects['y']):
		cond = (xmin < x < xmax) & (ymin < y < ymax)
		if cond:
			bool_aux = True
			# maskeo de estrella brillante
			cond2 = (178 < x < 190) & (189 < y < 201)
			if cond2:
				bool_aux = False
				break
	return bool_aux

def get_master(cube, filename):
	cube_array = fitsio.read(cube)
	master = np.median(cube_array, axis=0)
	hdu = fits.PrimaryHDU(master)
	hdulist = fits.HDUList([hdu])                # hdul = header data unit list
	hdulist.writeto('/home/consuelo/projects/FRBs/alopeke/reductions/' + filename + '.fits', overwrite = True)


def get_cube_reduced(cube, n_cube, master_bias, master_dark, flat):
	data = fitsio.read(cube)
	master_flat_bias = flat - master_bias
	#master_flat_dark = flat - master_dark

	reduced_bias = []
	#reduced_dark = []
	for i in range(1,5000):
		reduced_bias += [(data[i] - master_bias) / master_flat_bias]
		#reduced_dark += [(data[i] - master_dark) / master_flat_dark]

	reduced_bias = np.array(reduced_bias)
	#reduced_dark = np.array(reduced_dark)
	fits.writeto('/home/consuelo/projects/FRBs/alopeke/reductions/data_reduced/cube' + str(n_cube) + '_bias.fits', reduced_bias)
	#fits.writeto('/home/consuelo/projects/FRBs/alopeke/reductions/data_reduced/cube' + str(n_cube) + '_dark.fits', reduced_dark)


# for tables

def create_mask(x, y, sizebox):
	x = int(x)
	y = int(y)
	mask = np.zeros([256,256]).astype('int')
	mask[y-int(sizebox/2):y+int(sizebox/2), x-int(sizebox/2):x+int(sizebox/2)] = 1
	return mask

def bright_star_parameters(image, mask, rad_aperture_flux):
	image_masked = (image*mask).astype('int32')
	objects = sep.extract(image_masked, thresh=5, minarea=10)
	if len(objects) > 1:
		cond = np.where(objects['flux'] == np.max(objects['flux']))
		objects = objects[cond]
	flux, fluxerr, flag = sep.sum_circle(image_masked, objects['x'], objects['y'], rad_aperture_flux, gain=1.0)
	SN = flux/fluxerr
	return(flux, objects['x'], objects['y'], SN)
