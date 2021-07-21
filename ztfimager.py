#!/usr/bin/env python

"""
ztfimager queries ZTF and Pan-STARRS to plot both a 
ZTF reference image and multicolor Pan-STARRS PS1 image
to inspect visually inspect the field of view.
Author: 
    Joseph GUidry
For a description of updates, see the 
version_history.txt file.
Usage:
    ztfimage -ra|--ra -dec|--dec
Options:
    -ra --ra      RA of target in degrees
    -dec --dec    Declination of target in degrees
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

# Import local function
from image_funcs import ps_query, getztfrefurls, getps1colorim

# Suppress annoying warnings
import warnings
warnings.simplefilter('ignore')


#############################################################
## Generate arguments for command line parsing

parser = argparse.ArgumentParser()
parser.add_argument('-ra', '--ra',type=float,default=None,
                    help="Code name for telescope used.")
parser.add_argument('-dec', '--dec',type=float,default=None,
                    help="Code name for telescope used.")
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


#############################################################
## Perform PS1 query
ps1_query = ps_query(ra,dec,60)

## Get ZTF reference image
ztf_im = getztfrefurls(ra,dec)[0]
with fits.open(ztf_im) as hdul:
    ref_im = hdul[0].data
    hdr = hdul[0].header

## Get WCS coordinates from ZTF image header to plot locations of nearby stars found by Pan-STARRS
w = wcs.WCS(hdr)
target = [ra,dec]
coords1 = [ps1_query.raMean.values[0], ps1_query.decMean.values[0]]
coords2 = [ps1_query.raMean.values[1], ps1_query.decMean.values[1]]
coords3 = [ps1_query.raMean.values[2], ps1_query.decMean.values[2]]
wpix = w.wcs_world2pix(np.array([target,coords1,coords2,coords3]),1)

## Get the Pan-STARRS image and rotate to match ZTF orientation
ps1_im = np.flipud(rotate(ImageOps.flip(getps1colorim(ra,dec,size=500,filters="gri")),180))


#############################################################
## Plot the images
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6))

## Define nice scaling so the ZTF image is viewable:
ZS = ZScaleInterval(nsamples=10000, contrast=0.15, max_reject=0.5, 
                    min_npixels=5, krej=2.5, max_iterations=5)
vmin,vmax = ZS.get_limits(ref_im)

## Plot the images
ref_ax = ax1.imshow(np.flipud(ref_im), cmap='gray',origin='lower',vmin=vmin, vmax=vmax)
ax1.plot(30,30,'bx',ms=16,mew=4,alpha=0.7)
ax1.plot(wpix[3][0],wpix[3][1],'yx',ms=16,mew=4,alpha=0.4)
ax1.plot(wpix[2][0],wpix[2][1],'gx',ms=16,mew=4,alpha=0.4)
ax1.plot(wpix[1][0],wpix[1][1],'rx',ms=16,mew=4,alpha=0.9)
ax1.plot(wpix[0][0],wpix[0][1],'cx',ms=16,mew=4,alpha=0.6)
ax1.add_patch(Circle((wpix[0][0],wpix[0][1]), radius=7.0 ,
               edgecolor='c',facecolor='None',linewidth=1.75))
ax1.add_patch(Circle((30,30), radius=7.0 ,
               edgecolor='r',facecolor='None',linewidth=1.75))
ax1.set_title('ZTF Reference <filter> Image',fontsize=16)

ps1_ax =ax2.imshow(ps1_im, origin='lower')
ax2.set_title(r'PS1 Reference $gri$ Image',fontsize=16)
plt.show()
# plt.suptitle(wds[index],fontsize='xx-large',y=0.92);

# Display the PS1 query results

##########################################################
######### ADD PARSER ARG TO SHOW QUERY TABLE? ############
##########################################################


# ps1_query[['objID','raMean','decMean','nDetections','gMeanPSFMag','rMeanPSFMag','sep']]



# Print a closing message to the command line
print('\nFinished!\n')
