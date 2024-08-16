# Get image from the PS1 image server
# https://ps1images.stsci.edu/ps1image.html
# Add lines to download and write fits to save file.
# By Sophia Kim 2019.01.22. based on code PS1 suggests on the link above

# Pan-STARRS DR1 data query
# from https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/
# By CS Choi 


from __future__ import print_function
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from PIL import Image
from io import BytesIO
import numpy as np
import sys, os
import requests
import pylab

#size  = float(sys.argv[1]) # in pixel unit. 0.25"/pixel 2500 = 10', 3600 = 15', 5000=20', 6000=25'
#filt  = str(sys.argv[2]) 

#print('size = ',size)
#print('filter is ',filt)

# Helper functions to query the list of images and to extract images

def getimages(ra,dec,size=240,filters="grizy"):
    
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


def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="fits", color=False):
    
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
    table = getimages(ra,dec,size=size,filters=filters)
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


def getcolorim(ra, dec, size=240, output_size=None, filters="grizy", format="jpg"):
    
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
    url = geturl(ra,dec,size=size,filters=filters,output_size=output_size,format=format,color=True)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


def getgrayim(ra, dec, size=240, output_size=None, filter="g", format="fits"):
    
    """Get grayscale image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filter = string with filter to extract (one of grizy)
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    if filter not in list("grizy"):
        raise ValueError("filter must be one of grizy")
    url = geturl(ra,dec,size=size,filters=filter,output_size=output_size,format=format)
    r = requests.get(url[0])
    im = Image.open(BytesIO(r.content))
    return im



import requests 
from astropy.io.votable import parse_single_table 
 
def panstarrs_query(ra_deg, dec_deg, rad_deg, mindet=1, maxsources=10000, server=('https://archive.stsci.edu/'+'panstarrs/search.php')): 
    """
    Query Pan-STARRS DR1 @ MAST
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field 
                                          radius in degrees
                mindet: minimum number of detection (optional)
                maxsources: maximum number of sources
                server: servername
    returns: astropy.table object
    """
    r = requests.get(server, 
    params= {'RA': ra_deg, 'DEC': dec_deg, 
             'SR': rad_deg, 'max_records': maxsources, 
             'outputformat': 'VOTable', 
             'ndetections': ('>%d' % mindet)}) 
 
    # write query data into local file 
    outf = open('panstarrs.xml', 'w') 
    outf.write(r.text) 
    outf.close() 
 
    # parse local file into astropy.table object 
    data = parse_single_table('panstarrs.xml')
    return data.to_table(use_names_over_ids=True) 
 

#==========================================================================================
# To downlaod PS1 fits file, only use 'geturl' defined function.

#targetlist  = np.genfromtxt('/data2/sophia2/30inch/imsng_targetlist.txt', usecols=(0, 1, 2), dtype='str')
#targetlist = np.genfromtxt('target0426.list', usecols=(0, 1,2), dtype='str')
targetlist = np.genfromtxt('target0427.list', usecols=(0), dtype='str')
#save_dir    = '/data3/IMSNG/IMSNGgalaxies/refimg/PS1/'
#save_dir	  = '/data3/gwshare/refimg/PS1/'
save_dir 	= '/data7/sophia2/SN2023ixf/ref_frames/' #os.getcwd()
pixscale    = 0.4
tra, tdec = 210.9106542, +54.311675	 #2023ixf
#tra, tdec = 210.8024292, +54.34736111	#M101
#tra, tedc = 14.0533806, +54.3497500
size  = 3600
filt = "g"

#target_name = targetlist[:, 0]
#target_ra   = targetlist[:, 1]
#target_dec  = targetlist[:, 2]
fov	    = size*pixscale/60.
size = fov/pixscale*60.

filtlist = ['i', 'z', 'y']
for filt in filtlist:
	fitsurl=geturl(tra, tdec, size=int(size), filters=filt, format='fits')
	fh = fits.open(fitsurl[0])
	fh.writeto(f'./ref_frames/Ref-PS1-M101-{filt}-60arcmin.fits')

for file in targetlist:
	name=file.split('-')[2]
	data, header = fits.getdata(file, header=True)
	hdu_number=0 
	fits.getheader(file, hdu_number)
	ra1=header['RA'] 
	dec1=header['DEC']
	c    = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg))
	ra   = c.ra.degree
	dec  = c.dec.degree
	url  = geturl(ra, dec, size=size, output_size=None, filters=filt, format="fits")

	#if (save_dir+'Ref-PS1-/%s-%s-%darcmin.fits' %(name, filt, fov)) in: pass
	if os.path.isfile('Ref-PS1-%s-%s-%darcmin.fits' %(name, filt, fov))==True: pass
	else:
		try : 
			fh = fits.open(url[0])
			#fh.writeto(save_dir+'/%s/%s-%s-%darcmin-PS1.fits' %(filt, name, filt, fov))
			fh.writeto(save_dir+'/Ref-PS1-%s-%s-%darcmin.fits' %(name, filt, fov))
		except IndexError : error.append(i)	

print(target_name[error],' fail to download. \n Please check those images one by one')


for i in range(len(target_name)):
	name = target_name[i]
	ra1  = target_ra[i]
	dec1 = target_dec[i]
	c    = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg))
	ra   = c.ra.degree
	dec  = c.dec.degree
	url  = geturl(ra, dec, size=size, output_size=None, filters=filt, format="fits")

	#if (save_dir+'Ref-PS1-/%s-%s-%darcmin.fits' %(name, filt, fov)) in: pass
	if os.path.isfile('Ref-PS1-%s-%s-%darcmin.fits' %(name, filt, fov))==True: pass
	else:
		try : 
			fh = fits.open(url[0])
			#fh.writeto(save_dir+'/%s/%s-%s-%darcmin-PS1.fits' %(filt, name, filt, fov))
			fh.writeto(save_dir+'/Ref-PS1-%s-%s-%darcmin.fits' %(name, filt, fov))
		except IndexError : error.append(i)	


print(target_name[error],' fail to download. \n Please check those images one by one')

pixscale=0.25
filt	='r'
fov 	= size*pixscale/60.

error=[]
for i in range(len(data)):
	name	= data['host'][i]
	sn		= data['SN'][i]
	ra 		= data['ra'][i]
	dec		= data['dec'][i]
	url		= geturl(ra, dec, size=size, output_size=None, filters=filt, format='fits')

	try :
		fh	= fits.open(url[0])
		fh.writeto(save_dir+'/%s-host-%s-%s-%d.fits' %(sn, name, filt, fov))
	except IndexError : error.append(i)