#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import fits
import os
import csv
from astroquery.mast import Catalogs
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from astroquery.mast import Observations
from astropy.io import ascii

# In[18]:


def apass_query(ra_deg, 
                dec_deg, 
                rad_deg, 
                maxmag=18,
                minmag = 9,
                maxsources=100000):
    """
    Query APASS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    """
    vquery = Vizier(columns=['objID', 'RAJ2000', 'DEJ2000',
                             'e_RAJ2000', 'e_DEJ2000',
                             "Bmag", "e_Bmag",
                             "Vmag", "e_Vmag",
                             "g'mag", "e_g'mag",
                             "r'mag", "e_r'mag",
                             "i'mag", "e_i'mag"],
                    column_filters={"Bmag":
                                    ("<%f" % maxmag),
                                    "Vmag":
                                    ("<%f" % maxmag),
                                    "g'mag":
                                    ("<%f" % maxmag),
                                    "r'mag":
                                    ("<%f" % maxmag),
                                    "i'mag":
                                    ("<%f" % maxmag),
                                    "Bmag":
                                    (">%f" % minmag),
                                    "Vmag":
                                    (">%f" % minmag),
                                    "g'mag":
                                    (">%f" % minmag),
                                    "r'mag":
                                    (">%f" % minmag),
                                    "i'mag":
                                    (">%f" % minmag),

                                   },
                    
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='fk5')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="II/336/apass9")[0]


def sdss_query(ra_deg, dec_deg, rad_deg, maxmag=18,minmag = 9,maxsources=100000):
    """
    Query SDSS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    """
    vquery = Vizier(columns=['objID', 'RA_ICRS', 'DE_ICRS',
                             "umag", "e_umag",
                             "gmag", "e_gmag",
                             "rmag", "e_rmag",
                             "imag", "e_imag",
                             "zmag", "e_zmag",
                             "uc","gc","rc","ic","zc"],
                    column_filters={"umag":
                                    ("<%f" % maxmag),
                                    "gmag":
                                    ("<%f" % maxmag),
                                    "rmag":
                                    ("<%f" % maxmag),
                                    "imag":
                                    ("<%f" % maxmag),
                                    "zmag":
                                    ("<%f" % maxmag),
                                    "umag":
                                    (">%f" % minmag),
                                    "gmag":
                                    (">%f" % minmag),
                                    "rmag":
                                    (">%f" % minmag),
                                    "imag":
                                    (">%f" % minmag),
                                    "zmag":
                                    (">%f" % minmag)
                                   },
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="V/147/sdss12")[0]


def panstarrs1_query(ra_deg, dec_deg, rad_deg, maxmag=20,minmag = 13.5,
                    maxsources=1000000):
    """
    Query PanSTARRS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    """
    vquery = Vizier(columns=['*'],
                             #'RAJ2000', 'DEJ2000',
                             #'objID', 
                             #'e_RAJ2000', 'e_DEJ2000',
                             #'f_objID','Qual',
                             #'gmag', 'e_gmag',
                             #'rmag', 'e_rmag',
                             #'imag', 'e_imag',
                             #'zmag', 'e_zmag',
                             #'ymag', 'e_ymag',
                             #'gKmag', 'e_gKmag',
                             #'rKmag', 'e_rKmag',
                             #'iKmag', 'e_iKmag',
                             #'zKmag', 'e_zKmag',
                             #'yKmag', 'e_yKmag'],
                    column_filters={"gmag":
                                    ("<%f" % maxmag),
                                    "rmag":
                                    ("<%f" % maxmag),
                                    "imag":
                                    ("<%f" % maxmag),
                                    "gmag":
                                    (">%f" % minmag),
                                    "rmag":
                                    (">%f" % minmag),
                                    "imag":
                                    (">%f" % minmag),

                                   },
                    
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='fk5')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="II/349/ps1")[0]



def panstarrs2_query(ra,dec,radius):
    catalog_data = Catalogs.query_criteria(coordinates=f"{ra}{dec}", radius=radius,
                                            catalog="PANSTARRS", table="mean", data_release="dr2",
                                            nDetections=[("gte", 1)],
                                            qualityFlag=[('lte',64)],
                                            gMeanPSFMag=[('gte',13.5),('lte',21)],
                                            rMeanPSFMag=[('gte',13.5),('lte',21)],
                                            iMeanPSFMag=[('gte',13.5),('lte',21)],

                                            #zMeanPSFMag=[('gte',12),('lte',21)],
                                            #yMeanPSFMag=[('gte',12),('lte',21)],
                                            columns=["objID", "raMean","decMean","gMeanPSFMag","gMeanPSFMagErr","rMeanPSFMag","rMeanPSFMagErr","iMeanPSFMag","iMeanPSFMagErr","zMeanPSFMag","zMeanPSFMagErr","yMeanPSFMag","yMeanPSFMagErr"],
                                            sort_by=[("raMean")], pagesize=100000)
    return  catalog_data



def hms2dec(ra):
    RA = ra.split(':')
    h,m,s = float(RA[0]),float(RA[1]),float(RA[2])
    return 15*(h+m/60+s/3600)
def dms2dec(dec):
    Dec = dec.split(':')
    if dec[0]== '-':
        d,m,s = float(Dec[0]),float(Dec[1]),float(Dec[2])
        return (d-m/60-s/3600)
    elif dec[0]== '+':
        d,m,s = float(Dec[0]),float(Dec[1]),float(Dec[2])
        return (d+m/60+s/3600)





#%%
from astropy.table import Table
Targetlist = Table()
Targetlist['obj'] = ['NGC1784']
Targetlist['ra'] = ['05:05:27']
Targetlist['dec'] = ['-11:52:17']
FIELD = Targetlist['obj']
RA = Targetlist['ra']
DEC = Targetlist['dec']

Field = []
for i in range(len(FIELD)):
    ra = hms2dec(RA[i])
    dec=dms2dec(DEC[i])
    Field.append([FIELD[i],ra,dec])

#%%
######### CHANGE DIRECTORY HERE ##########
#Field = os.listdir('/data2/hhchoi1022/IMSNG/SOAO/Data')
#os.chdir('/data2/hhchoi1022/IMSNG/SOAO/SkyCatalog/')

Targetlist = ascii.read('/home/hhchoi1022/Desktop/Gitrepo/config/alltarget.dat')
FIELD = Targetlist['obj']
RA = Targetlist['ra']
DEC = Targetlist['dec']

Field = []
for i in range(len(FIELD)):
    ra = hms2dec(RA[i])
    dec=dms2dec(DEC[i])
    Field.append([FIELD[i],ra,dec])
#%%

os.chdir('/data1/Skycatalog/PanSTARRS1')
for field,ra,dec in Field:
    try:
        tab = panstarrs1_query(ra,dec,0.5)
        table = tab.to_pandas()
        if len(table)>10:
            table.to_csv(f'{field}.csv',
                    sep = ',',
                    na_rep ='',
                    index =False)
            print(f'{field} Downloaded')
        else:
            print(f'{field} not covered by Panstarrs')
    except:
        print(f'SOMETHING WRONG for {field}')
        




# %%
