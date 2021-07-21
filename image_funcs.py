'''
Script originally written by Zach Vanderbosch (https://github.com/zvanderbosch)

Updated by Joseph Guidry to sort outputted DataFrame by proximity to the queried target.
'''

from __future__ import print_function
import io
import os
import sys
import copy
import json
import time
import re
import requests
import warnings
from io import StringIO, BytesIO
import numpy as np
import pandas as pd
from glob import glob
from astropy.io import fits,ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.time import Time
from astropy.table import Table
import pylab
from PIL import Image

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib 



def ps_query(ra,dec,rad,full_query=False):

    # Create an Empty DataFrame
    # if full_query:
    columns = """objID,raMean,decMean,nDetections,ng,nr,ni,nz,ny,
            gMeanPSFMag,gMeanPSFMagErr,rMeanPSFMag,rMeanPSFMagErr,
            iMeanPSFMag,iMeanPSFMagErr,zMeanPSFMag,zMeanPSFMagErr,
            yMeanPSFMag,yMeanPSFMagErr,gMeanKronMag,gMeanKronMagErr,
            rMeanKronMag,rMeanKronMagErr,iMeanKronMag,iMeanKronMagErr,
            zMeanKronMag,zMeanKronMagErr,yMeanKronMag,yMeanKronMagErr""".split(',')
    # else:
    #     columns = """objID,raMean,decMean,gMeanPSFMag,rMeanPSFMag""".split(',')
    columns = [x.strip() for x in columns]
    columns = [x for x in columns if x and not x.startswith('#')]

    # Query Constraints
    qpath = './'
    constraints = {'ng.gt':4,  # Minimum 5 Detections per filter
                   'nr.gt':4,
                   'ni.gt':4,
                   'nz.gt':4,
                   'ny.gt':4}
    # constraints = {'nDetections.gt':5}

    # Define central RA/Dec of the image
    ra_center = ra
    dec_center = dec
    # Search Radius Defined to Reach From Center-to-Corner of Image
    radius = rad/3600.0

    # Perform the Cone Search
    results = ps1cone(ra_center,dec_center,radius,release='dr2',columns=columns,**constraints)

    # Convert Results into an Astropy Table, improve formatting,
    # and then create a Pandas Table
    apy_tab = ascii.read(results)
    for f in 'grizy':
    # for f in 'gr':
        col1 = f+'MeanPSFMag'
        col2 = f+'MeanPSFMagErr'
        try:
            apy_tab[col1].format = ".4f"
            apy_tab[col1][apy_tab[col1] == -999.0] = np.nan
        except KeyError:
            print("{} not found".format(col1))
        try:
            apy_tab[col2].format = ".4f"
            apy_tab[col2][apy_tab[col2] == -999.0] = np.nan
        except KeyError:
            print("{} not found".format(col2))
    ps_tab = apy_tab.to_pandas()

    # Sort the DataFrame to ensure the target object occupies the first row
    sep = []
    c0 = SkyCoord(ra*u.deg,dec*u.degree,frame='icrs')
    for r,d in zip(ps_tab.raMean.values,ps_tab.decMean.values):
        c = SkyCoord(r*u.deg,d*u.degree,frame='icrs')
        sep.append(c0.separation(c).arcsecond)
    sep = np.array(sep)
    sep -= np.sort(sep)[0] 
    num_cols = len(ps_tab.columns)
    ps_tab.insert(num_cols+0,'sep',sep)
    sorted_tab = ps_tab.sort_values('sep',ascending=True).reset_index(drop=True) 
    
    # Save the query for future use
    # ps_tab.to_csv(qpath+'ps_query.csv',index=False)

    return sorted_tab



"""
This script contains several functions useful for
retrieving a variety of ZTF data products such as:
    - Reference images (from IPAC)
    - Science images (from IPAC)
    - Light Curves (from IPAC)
    - Transient Alert Packets (from MARS/LCO)

Author: Zach Vanderbosch
"""


# Function to retrieve URLs for reference image cutouts
# centered on a given RA and DEC
def getztfrefurls(ra,dec,size=45):
    """
        ra   = Right Ascension in Decimal Degrees
        dec  = Declination in Decimal Degrees
        size = Image size in arcseonds

    Returns:
        urls = list of URLs used to download image data
    """
    
    # First get some info related to the reference image
    search_url = 'https://irsa.ipac.caltech.edu/ibe/search/ztf/products/ref'
    pos_url = search_url + '?POS={:.3f},{:.3f}'.format(ra,dec)
    r1 = requests.get(pos_url)
    im_table = ascii.read(r1.text).to_pandas()

    # Get image meta data
    urls = []
    maxbit = 33554432 # Quality cut on the INFOBITS parameter
    num_image = len(im_table[im_table.infobits < maxbit]) # Number of images in table
    for im in range(num_image):
        field = str(im_table[im_table.infobits < maxbit].field.iloc[im]).zfill(6)
        filtcode = str(im_table[im_table.infobits < maxbit].filtercode.iloc[im])
        ccdid = str(im_table[im_table.infobits < maxbit].ccdid.iloc[im]).zfill(2)
        qid = str(im_table[im_table.infobits < maxbit].qid.iloc[im])

        # Construct the Image Download URL
        data_url = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/ref'
        spec_url = '/{}/field{}/{}/ccd{}/q{}/'.format(field[0:3],field,filtcode,ccdid,qid)
        imname = 'ztf_{}_{}_c{}_q{}_refimg.fits'.format(field,filtcode,ccdid,qid)
        condis = '?center={:.5f},{:.5f}&size={}arcsec&gzip=false'.format(ra,dec,size)
        urls.append(data_url + spec_url + imname + condis)

    return urls


'''
Functions written by MAST as a part of the Pan-STARRS API
for useful retrieval.

Sourced from: http://ps1images.stsci.edu/ps1_dr2_api.html
'''




# Function to calculate angular separation between two RA/Dec coordinates
def angle_sep(d1,d2,a1,a2):
    rad = np.radians(np.asarray([d1,d2,a1,a2]))
    asep = np.arccos(np.sin(rad[0])*np.sin(rad[1])+
                     np.cos(rad[0])*np.cos(rad[1])*
                     np.cos(rad[3]-rad[2]))
    return np.degrees(asep)*3600.0 # Separation in Arcseconds 


def ps1search(table="mean",release="dr1",format="csv",columns=None,
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
           **kw):
    """Do a general search of the PS1 catalog (possibly without ra/dec/radius)
    
    Parameters
    ----------
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'nDetections.min':2).  Note this is required!
    """
    
    data = kw.copy()
    if not data:
        raise ValueError("You must specify some parameters for search")
    checklegal(table,release)
    if format not in ("csv","votable","json"):
        raise ValueError("Bad value for format")
    url = "{baseurl}/{release}/{table}.{format}".format(**locals())
    if columns:
        # check that column values are legal
        # create a dictionary to speed this up
        dcols = {}
        for col in ps1metadata(table,release)['name']:
            dcols[col.lower()] = 1
        badcols = []
        for col in columns:
            if col.lower().strip() not in dcols:
                badcols.append(col)
        if badcols:
            raise ValueError('Some columns not found in table: {}'.format(', '.join(badcols)))
        # two different ways to specify a list of column values in the API
        # data['columns'] = columns
        data['columns'] = '[{}]'.format(','.join(columns))

# either get or post works
#    r = requests.post(url, data=data)
    r = requests.get(url, params=data, timeout=500)

    if verbose:
        print(r.url)
    r.raise_for_status()
    if format == "json":
        return r.json()
    else:
        return r.text

def checklegal(table,release):
    """Checks if this combination of table and release is acceptable
    
    Raises a VelueError exception if there is problem
    """
    
    releaselist = ("dr1", "dr2")
    if release not in ("dr1","dr2"):
        raise ValueError("Bad value for release (must be one of {})".format(', '.join(releaselist)))
    if release=="dr1":
        tablelist = ("mean", "stack")
    else:
        tablelist = ("mean", "stack", "detection")
    if table not in tablelist:
        raise ValueError("Bad value for table (for {} must be one of {})".format(release, ", ".join(tablelist)))


def ps1metadata(table="mean",release="dr1",
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"):
    """Return metadata for the specified catalog and table
    
    Parameters
    ----------
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    baseurl: base URL for the request
    
    Returns an astropy table with columns name, type, description
    """
    
    checklegal(table,release)
    url = "{baseurl}/{release}/{table}/metadata".format(**locals())
    r = requests.get(url)
    r.raise_for_status()
    v = r.json()
    # convert to astropy table
    tab = Table(rows=[(x['name'],x['type'],x['description']) for x in v],
               names=('name','type','description'))
    return tab


def ps1cone(ra,dec,radius,table="mean",release="dr1",format="csv",columns=None,
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
           **kw):
    """Do a cone search of the PS1 catalog
    
    Parameters
    ----------
    ra (float): (degrees) J2000 Right Ascension
    dec (float): (degrees) J2000 Declination
    radius (float): (degrees) Search radius (<= 0.5 degrees)
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'nDetections.min':2)
    """
    
    data = kw.copy()
    data['ra'] = ra
    data['dec'] = dec
    data['radius'] = radius
    return ps1search(table=table,release=release,format=format,columns=columns,
                    baseurl=baseurl, verbose=verbose, **data)

###### Additional Functions sourced from https://ps1images.stsci.edu/ps1image.html ########

# Query Pan-STARRS
def getps1images(ra,dec,size=240,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table

def getps1url(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getps1images(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url

def getps1colorim(ra, dec, size=240, output_size=None, filters="grizy", format="jpg"):
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    url = getps1url(ra,dec,size=size,filters=filters,output_size=output_size,format=format,color=True)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


