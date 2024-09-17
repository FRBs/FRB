###### README #####

# author: Clancy W. James
# email: clancy.james@curtin.edu.au

# Purpose: perform a single statistical test:
# does fitting a mean reduce the variance? But it does this using weights according to estimated errors

# Run using:
# python weighted_multi_image_fit.py.py infile1 infile2 ...
# e.g. python multi_image_fit.py 190102.dat 190608.dat 190711.dat
#     Infiles have: one column per source, rows are ra_offset, dec_offset, ra_err, dec_err
#     Please use simple text file with numbers and whitespace separation

# Output: some statistical shit

import numpy as np
import sys
import os
import scipy.stats as stats
from matplotlib import pyplot as plt


def load_data(infile):
	print("Loading data infile ",infile,"...")
	if os.path.isfile(infile):
		try:
			data=np.loadtxt(infile) # atual data read statement
		except:
			print("ERROR: error reading infile ",infile)
			exit()
	else:
		print("ERROR: Could not locate input file ",infile)
		exit()
	dims=data.shape
	
	if len(dims) != 2 or dims[0] != 4:
		print("ERROR: Please ensure data is in correct format")
		print("Four rows: ra offset, dec offset, ra err, dec err")
		print("and one column per source.")
		exit()
	print("      Success! Found ",dims[1]," sources")
	return data


# Main program, one big module
def main():
	
	print("\n\n\n\n")
	print("          #########################################################")
	print("          ## Clancy's slightly less dodgy error analysis program ##")
	print("          #########################################################\n\n")
	
	
	
	########### Checking all is OK #########
	nargs=len(sys.argv)
	if nargs == 1:
		print("ERROR: Please enter input file type as first argument")
		exit()
	
	nfiles=nargs-1
	
	print("Reading in data from ",nfiles," input files")
	data=[]
	ndtot=0
	ndoff=0
	for i in np.arange(nfiles):
		data.append(load_data(sys.argv[i+1]))
		
	
	print("\n\n\n\n ########## H0: no offsets ########")
	# This part of the routine tests the H0 that the measured positions and errors
	# are consistent with 0,0 at the stated level of uncertainty
	# it simply does this by calculating the variance in the offsets, and in
	# the stated errors.
	
	### estimates global sigmas ###
	# estimated variance based on offsets
	var=0
	var_dec=0
	var_ra=0
	nsrc=np.zeros([nfiles])
	var_w=0
	
	######## loops through sources #######
	# artificial scale
	for i,frb in enumerate(data):
		nsrc[i]=(frb.shape[1])
		# weight each point with 1/err**2
		# we weight the  measured offsets by 1/err.
		# Under H0, this makes them come from a standard normal distribution
		wr=frb[0,:]/frb[2,:]
		wd=frb[1,:]/frb[3,:]
		
		# we now calculate the variances using the std normal deviations
		vr = np.sum(wr**2)
		vd = np.sum(wd**2)
		# add these to cumulative sums over all the FRB sets
		
		var += vr+vd
		var_ra += vr
		var_dec += vd
		
		# I do not know what this is for
		var_w += np.sum((1./frb[2,:])**2)+np.sum((1./frb[3,:])**2)
		
		
		n=frb[0,:].size
		meanr=np.sum(frb[0,:])/n
		meand=np.sum(frb[1,:])/n
		offr=np.abs(frb[0,:]-meanr)
		offd=np.abs(frb[1,:]-meand)
		
		plt.figure()
		plt.xlabel('Offet from unweighted mean')
		plt.ylabel('Error on point')
		plt.plot(offr,frb[2,:],color='red',linestyle='',marker='x',label='ra')
		plt.plot(offd,frb[3,:],color='blue',linestyle='',marker='o',label='dec')
		plt.legend()
		plt.tight_layout()
		plt.savefig('err_vs_offset_'+str(i)+'.pdf')
		plt.close()
	
	# Under H0, var, var_ra, and var_dec should come from a chi2 distribution with
	# nsrc_tot degrees of freedom
	
	# calculate the degrees of freedom for var_ra (nsrc) and var (2*this due to ra and dec)
	nsrc_tot=np.sum(nsrc)
	ndat_tot=nsrc_tot*2
	
	# Under H0, the following should be ~1. This helps estimate how much larger/smaller
	# the true number must be
	sig_ra=(var_ra/nsrc_tot)**0.5
	sig_dec=(var_dec/nsrc_tot)**0.5
	sig=(var/nsrc_tot/2.)**0.5
	
	############ Testing that variance between ra and dec is equal IF H0 is true ###########
	# We now perform an f-test for equality of variance between ra and dec
	# the goal is to work out what the errors would have to be in oder to satisfy H0
	f_stat=var_ra/var_dec
	#NOTE: use nsrc_tot, not nsrc_tot-1, since we *know* the means in this example - it is 0,0 (H0)
	cdf=stats.f.cdf(f_stat,nsrc_tot,nsrc_tot)
	# performs a two-sided test for the rtio being too small or too big.
	p_val=min(cdf,1.-cdf)*2.
	print("Under H0, errors in ra,dec are {0:6.2f} {1:6.2f} greater than quoted in file".format(sig_ra,sig_dec))
	print("Assuming they are intrinsically equal (two-sided p-value of {0:6.2f}), this factor is {1:6.2f}".format(p_val,sig))
	
	
	print("\n\n\n\n ########## H1: offsets ########")
	### now calculates mean offsets for each source ###
	
	mean_ra=np.zeros([nfiles])
	mean_dec=np.zeros([nfiles])
	mean_ra_err=np.zeros([nfiles])
	mean_dec_err=np.zeros([nfiles])
	var_ra_mean=0
	var_dec_mean=0
	var_mean=0
	print("Calculating offsets (assuming correct errors)...")
	print("Source  ra_off (u/w scat err) [w scat err] [[input error]]     dec_off (u/w scat err) [w scat err] [[input error]]:")
	for i,frb in enumerate(data):
		
		# This returns an estimate for the weighted mean and an error in that mean
		# i.e. it is reduced by the sqrt N factor
		# previously this was bugged, and did not use the optimal weighting (1/err^2)
		mean_ra[i],mean_ra_err[i]=get_weighted_mean(frb[0,:],frb[2,:])
		mean_dec[i],mean_dec_err[i]=get_weighted_mean(frb[1,:],frb[3,:])
		
		# number of degrees of freedom in one dimension
		npts=frb[0,:].size
		ndf=npts-1. #for each coordinate
		reldf=ndf/npts
		
		#### standard calculations ######
		# we now calculate the offsets given our estimate for the mean offset
		dwr=(frb[0,:]-mean_ra[i])
		dwd=(frb[1,:]-mean_dec[i])
		
		# calculates unweighted scatter error
		# first lines cal the std dev of the scatter
		# next the error in the mean
		sigr=(np.sum(dwr**2)/ndf)**0.5
		sigd=(np.sum(dwd**2)/ndf)**0.5
		sigr /= npts**0.5
		sigd /= npts*0.5
		
		# if the errors are correct, then the below should be a standard normal deviate
		# i.e. any deviation from 0 is purely due to random variation
		wr=dwr/frb[2,:]
		wd=dwd/frb[3,:]
		
		###### calculations for chi2 tests only #######
		# calculates new variance under H0: no true offset
		# will come from a chi2 distribution with (2*) n-nfrb degrees of freedom
		vr = np.sum(wr**2)
		vd = np.sum(wd**2)
		# we do this for ra and dec separately, and together
		var_mean += vr+vd
		var_ra_mean += vr
		var_dec_mean += vd
		
		######### Calculations for error in mean ########
		# reduces this for the weighting, i.e. so it's agnostic to the values of the weights
		# this is identical to the std dev of the weighted mean, but re-normalised by the
		# sum of the weights, such that it is only the *relative* weights that count
		# the base offsets are also derived from data, not the input values of std
		# total formula is \sqrt{ (\sum_i w (xi-xbar)^2 )/((M-1)/M * \sum_i w_i) }
		# here the numerator is simply vr and vd
		# this gives the standard deviation of the sample. The error in the mean 
		# will be sqrt(n) less than this, i.e. the factor M-1/M becomes M-1
		wvr = vr/(ndf*np.sum(frb[2,:]**-2))
		wvd = vd/(ndf*np.sum(frb[3,:]**-2))
		
		errr = wvr**0.5
		errd = wvd**0.5
		print("   {0:d} {1:8.2f}    ({2:6.2f})       [{3:6.2f}]   [[{4:6.2f}]]  {5:14.2f}    ({6:6.2f})      [{7:6.2f}]    [[{8:6.2f}]]".format(i,mean_ra[i],sigr,errr,mean_ra_err[i],mean_dec[i],sigd,errd,mean_dec_err[i]))
		#mean_dec[i],mean_ra_err[i],mean_dec_err[i],errr,errd))
		#print("  The chi2 from subsequent fits is ",errr,errd)
	
	print("In the above, 'scatter error' is derived by ignoring the quoted values of")
	print("the errors, and simply looking at the scatter of points about the mean.")
	print("Input error is estimated by assuming that the quoted errors in the input file are correct.")
	print("The scatter error is more robust.")
	print("The input error is more rigourously correct IF the inputs are sensible.")
	print("In general, take the input error unless the scatter error is significantly different.")
	
	all_means=np.concatenate((mean_ra,mean_dec))
	all_errs=np.concatenate((mean_ra_err,mean_dec_err))
	sigma_offsets=(np.sum(all_means**2/all_errs)/np.sum(1./all_errs))**0.5
	#print("Standard deviation of offsets about 0 is {0:6.2f}".format(sigma_offsets))
	
	# performs an f-test for reduction in variance between 0,0 model and offsets in ra,dec
	# NOTE: if setting ra OR dec offset to zero, then the result will be different
	print("\n\n\n\n#### Testing for rejection of H0 ####")
	# test for reduction in variance
	# calculates f-stat for Chow test
	ndf1=nsrc_tot*2 # degrees of freedom for total number of points
	ndf2=2*nsrc_tot-2*nfiles # df including means
	f_stat=( (var-var_mean)/ndf1 )  / ( var_mean / ndf2 )
	cdf=stats.f.cdf(f_stat,ndf1,ndf2)
	# reject if F is large
	print("Fitting offsets reduced variance by a factor of {0:6.2f}".format(var/var_mean))
	print("Probability under the no-offset hypothesis that fitting")
	print("     offsets would lead to at least this much reduction ")
	print("     in variance is {0:5.2f}".format(1.-cdf))
	
	
	print("\n\n\n\n ##### Is the new offset consistent with the errors? \n")
	ndf_1d=ndf2/2 # for ra or dec only
	p_both=1.-stats.chi2.cdf(var_mean,ndf2)
	p_ra=1.-stats.chi2.cdf(var_ra_mean,ndf_1d)
	p_dec=1.-stats.chi2.cdf(var_dec_mean,ndf_1d)
	print("Coord     Variance     ndf       p \n")
	print("Both: {0:10.2f} {1:4d}  {2:7.2E}".format(var_mean,int(ndf2),p_both))
	print("  Ra: {0:10.2f} {1:4d}  {2:7.2E}".format(var_ra_mean,int(ndf_1d),p_ra))
	print(" Dec: {0:10.2f} {1:4d}  {2:7.2E}".format(var_dec_mean,int(ndf_1d),p_dec))
	print("(If p-values are low, then random errors will be larger than given in inputs.)")
	print("(If p-values are high, then random errors are consistent with inputs.)")
	
	
	
	# we now estimate the error on the mean using observed differences
	print("\n\n\n ##### Estimating new mean from observed offsets ####\n")
	
	crit_90=stats.f.ppf(0.9,ndf1,ndf2)
	var_mean_90=(var/ndf1) / (1./ndf1 + crit_90/ndf2)
	std_dev_offsets=((var-var_mean_90)/var_w)**0.5
	
	print("\n\n\n\b### Estimating characteristic offset scale which could have been detected #####")
	print("At a 90% confidence level, the critical reduction in variance would be {0:6.2f}".format(var/var_mean_90))
	print("   This corresponds to an average std deviation in offsets of {0:8.2f}".format(std_dev_offsets))
	print("   CAUTION: this is only an estimate, not a hard statistical statement")
	
# calculates weighted means and error on them from list of observations and errors
def get_weighted_mean(obs,err):
	# we use the weight as 1/err**2
	mean=np.sum(obs/err**2)/np.sum(1./err**2)
	# the std dev in this mean is 
	err=(np.sum(1./err**2))**-0.5
	#err=obs.size**0.5/np.sum(1./err**2)
	return mean,err
	

	
main()
