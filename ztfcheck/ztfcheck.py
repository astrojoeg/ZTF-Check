#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ztfimager queries ZTF and Pan-STARRS to plot both a 
ZTF reference image and multicolor Pan-STARRS PS1 image
to inspect visually inspect the field of view.
Author: 
    Joseph GUidry
For a description of updates, see the 
version_history.txt file.
Usage:
    ztfimage -ra|--ra -dec|--dec -q|--query
Options:
    -ra --ra      RA of target in degrees
    -dec --dec    Declination of target in degrees
    -q --query	  Print out the table of results for neighboring stars
    			  within the Pan-STARRS query
""" 


# Import necessary packages
import argparse
from astropy import wcs
from astropy.io import fits
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
from PIL import Image, ImageOps
from scipy.ndimage import rotate
import sys

# Import local functions
from image_funcs import ps_query, getztfrefurls, getps1colorim

# Suppress annoying warnings
import warnings
warnings.simplefilter('ignore')


#############################################################
## Define "main" function to make ZTF Check a command line executable
def main():
	## Generate arguments for command line parsing
	parser = argparse.ArgumentParser()
	parser.add_argument('-ra', '--ra',type=float,default=None,
	                    help="Right ascension of target.")
	parser.add_argument('-dec', '--dec',type=float,default=None,
	                    help="Declination of target.")
	parser.add_argument('-q', '--query',action='store_true',
	                    help="Whether to print the results from the 30 arcsec Pan-STARRS query.")
	args = parser.parse_args()
	ra = args.ra 
	dec = args.dec 

	## Check that inputted ra and dec are valid:
	if ra == None or type(ra) != float or ra < 0. or ra > 360.:
		print('\nERROR! "{}" is not a valid RA input.'.format(ra))
		print('Supported RA values in degrees are:')
		print('0. < RA < 360.')
		print("")
		sys.exit(1)

	if dec == None or type(dec) != float or dec < -90. or dec > 90.:
		print('\nERROR! "{}" is not a valid decl. input.'.format(dec))
		print('Supported decl. values in degrees are:')
		print('-90. < RA < 90.')
		print("")
		sys.exit(1)

	print("")	


	#############################################################
	## Perform PS1 query
	try:
		ps1_query = ps_query(ra,dec,30,full_query=False)
		if args.query:
			print("Results from 30 arcsec Pan-STARRS query:\n")
			print(ps1_query[['raMean','decMean','gMeanPSFMag','rMeanPSFMag','separation']])
			print("")
	except FileNotFoundError:
		print('ERROR! "R.A. = {}, decl. = {}" is not covered by Pan-STARRS.'.format(ra, dec))
		print('Please check your coordinates or try another target.')
		print("")
		sys.exit(1)

	## Get ZTF reference image
	try:
		ztf_im = getztfrefurls(ra,dec)[0]
		filt = ztf_im.split('/')[10]
		with fits.open(ztf_im) as hdul:
		    ref_im = hdul[0].data
		    hdr = hdul[0].header
		    # print(hdr)
	except IndexError:
		print('ERROR! "R.A. = {}, decl. = {}" is not covered by ZTF.'.format(ra, dec))
		print('Please check your coordinates or try another target.')
		print("")
		sys.exit(1)

	## Get WCS coordinates from ZTF image header to plot locations of nearby stars found by Pan-STARRS
	w = wcs.WCS(hdr)
	target = [ra,dec]
	neighbor1 = [ps1_query.raMean.values[1], ps1_query.decMean.values[1]]
	neighbor2 = [ps1_query.raMean.values[2], ps1_query.decMean.values[2]]
	neighbor3 = [ps1_query.raMean.values[3], ps1_query.decMean.values[3]]
	wpix = w.wcs_world2pix(np.array([target,neighbor1,neighbor2,neighbor3]),1)

	## Get the Pan-STARRS image and rotate to match ZTF orientation
	ps1_im = np.flipud(rotate(ImageOps.flip(getps1colorim(ra,dec,size=500,filters="gri")),180))


	#############################################################
	## Plot the images
	fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6), tight_layout=True, gridspec_kw={'wspace':0.25})

	## Define nice scaling so the ZTF image is viewable:
	ZS = ZScaleInterval(nsamples=10000, contrast=0.15, max_reject=0.5, 
	                    min_npixels=5, krej=2.5, max_iterations=5)
	vmin,vmax = ZS.get_limits(ref_im)

	## Plot the images
	ref_ax = ax1.imshow(np.flipud(ref_im), cmap='gray',origin='lower',vmin=vmin, vmax=vmax)
	ax1.plot(wpix[3][0],wpix[3][1],'C1x',ms=16,mew=4,alpha=1,label='Neighbor 3')
	ax1.plot(wpix[2][0],wpix[2][1],'C2x',ms=16,mew=4,alpha=1,label='Neighbor 2')
	ax1.plot(wpix[1][0],wpix[1][1],'C3x',ms=16,mew=4,alpha=1,label='Neighbor 1')
	ax1.plot(wpix[0][0],wpix[0][1],'C9x',ms=16,mew=4,alpha=1,label='Target')
	ax1.add_patch(Circle((wpix[0][0],wpix[0][1]), radius=5.0 ,
	               edgecolor='C9',facecolor='None',linewidth=1.75,label='5 arcsec Aperture'))
	ax1.set_title(r'ZTF Reference ${}$-band Image'.format(filt[1]),fontsize=16)
	ax1.tick_params(which='major',direction='out',labelsize=12,width=1.25,length=6)
	ax1.set_xlabel('Separation (arcsec)',fontsize=14)
	ax1.set_ylabel('Separation (arcsec)',fontsize=14)
	ax1.legend(loc='best',markerscale=0.70)

	ps1_ax =ax2.imshow(ps1_im, origin='lower')
	ax2.set_title(r'PS1 Reference $gri$ Image',fontsize=16)
	ax2.tick_params(which='major',direction='out',labelsize=12,width=1.25,length=6)
	ax2.set_xlabel('Separation (arcsec)',fontsize=14)
	ax2.set_ylabel('Separation (arcsec)',fontsize=14)

	plt.show()

	## Print a closing message to the command line
	print('\nFinished!\n')
