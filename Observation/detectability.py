#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 22:08:40 2022

@author: hhchoi1022
"""
conf_param = dict(
                  obspath          = "/home/hhchoi1022/Desktop/Gitrepo/config/obs_location.txt",
                  catpath          = '/home/hhchoi1022/Desktop/Gitrepo/config/alltarget.dat',
                  savepath         = '/home/hhchoi1022/Desktop/Gitrepo/observatory/KCT/ACP/',
                  observatory      = 'KCT',
                  altlimit         = 30,
                  sunlimit         = -18,
                  moonseperation   = 40,
                  moonphase = True,
                  frame = 'fk5',
                  graph = True
                  )

#%%
import gcn
import lxml.etree
from astropy.io import ascii
import numpy as np
import pytz
import ephem
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
from astropy.visualization import astropy_mpl_style, quantity_support
import matplotlib.pyplot as plt
from rtsmaker import rtsmaker
plt.style.use(astropy_mpl_style)
quantity_support()


def parse_xml(xmlfile):
    payload = open(xmlfile, 'rb').read()
    root = lxml.etree.fromstring(payload)
    pos2d = root.find('.//{*}Position2D')
    when = root.find(".//{*}Time")
    ##### Set required parameters for output #####
    output = dict(
                  ID = str(root.find(".//Param[@name='TrigID']").attrib['value']),
                  ra = float(pos2d.find('.//{*}C1').text),
                  dec = float(pos2d.find('.//{*}C2').text),
                  isotime = when.find(".//{*}ISOTime").text
                  )
    return output 

def check_detectability(ra, dec, observatory, time, graph = True, frame = 'fk5', sunlimit = -18, altlimit = 30, ID = None):
    # input  : ra, dec, observatory information, time(UT)
    # output : detectability graph, obsstarttime, booltype(detectability)
    # Source code and graph design by Gregory. S.H. Paek
    # Modified by Hyeonho Choi
    
    # Observatory information
    obsinfo          = ascii.read(conf_param['obspath'])
    obsname          = np.copy(obsinfo['name'])
    obsindex         = np.where(obsname == observatory)[0]
    obslat           = (np.copy(obsinfo['latitude(N+)'])[obsindex])[0]
    obslon           = (np.copy(obsinfo['longitude(E+)'])[obsindex])[0]
    obsalt           = (np.copy(obsinfo['altitude'])[obsindex])[0]
    loc = EarthLocation(lat=obslat*u.deg, lon=obslon*u.deg, height=obsalt*u.m)
    
    # Time information
    t = Time(time, format='isot', scale='utc')
    delta_midnight = np.linspace(0, 24, 200)*u.hour
    times = t + delta_midnight
    earthframe = AltAz(obstime=times, location=loc)
    
    # Query the object
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame=frame)
    c_altaz = c.transform_to(earthframe)
    
    # Query the Sun
    c_sun_altaz = get_sun(times).transform_to(earthframe)
    indx_sun = np.where(
        (c_sun_altaz.alt < sunlimit*u.deg)
        )
    indx_obj = np.where(
        (c_sun_altaz.alt < sunlimit*u.deg)&
        (c_altaz.alt > altlimit*u.deg)
        )
    
    # Query the Moon
    moonseperation = conf_param['moonseperation']
    c_moon_altaz = get_moon(times).transform_to(earthframe)
    indx_moon, sep_moon, _ = c_altaz.match_to_catalog_sky(c_moon_altaz) #  Moon distance
    sep_moon_sun = sep_moon[indx_sun][sep_moon[indx_sun]>moonseperation*u.deg]  #   Moon distance (@night)
    moondist_tar = np.median(sep_moon_sun)
    mphase       = ephem.Moon(t.iso)
   
    # Object detectability
    if len(indx_obj[0])==0:
        fir_time = times[0]
        fir_alt = 0
        key_detect = 0
    else:
        fir_time = np.min(times[indx_obj])
        fir_alt = c_altaz.alt[indx_obj[0][np.argmin(times[indx_obj])]]
        key_detect = (fir_time-t).jd*24<12
        #key_detect = (fir_time-t).jd*24<12
        if not key_detect:
            key_detect = 2 # For the case that the object is detectable at next night 
    
    # Visualization
    if graph:
        plt.close('all')
        plt.rcParams["figure.figsize"] = (15,6)
        plt.plot(delta_midnight, c_sun_altaz.alt, color='tomato', label='Sun', alpha=0.5)
        plt.plot(delta_midnight, c_moon_altaz.alt, color='dodgerblue', ls='--', label='Moon', alpha=0.5)
        plt.scatter(delta_midnight, c_altaz.alt, c=c_altaz.az, label='Target', marker = 'o', lw=0, s=25, cmap='viridis')
        clb = plt.colorbar()
        clb.set_label('Azimuth [deg]')
        plt.fill_between(delta_midnight, -90*u.deg, 90*u.deg, c_sun_altaz.alt < -0*u.deg, color='0.75', zorder=0)
        plt.fill_between(delta_midnight, -90*u.deg, 90*u.deg, c_sun_altaz.alt < -18*u.deg, color='0.65', zorder=0)
        plt.text(1, 85, f'Moon seperation = {int(moondist_tar.degree)}', fontsize = 15)
        plt.text(1, 80, 'Moon phase = %.2f' %mphase.moon_phase, fontsize = 15)
        plt.legend(loc='upper center', fontsize=16, framealpha=1.0)
        plt.xlim([np.min(delta_midnight), np.max(delta_midnight)])
        plt.xticks(np.arange(0, 24+1))
        # plt.xlabel('Hours from EDT Midnight', fontsize=20)
        plt.xlabel(f'Hours after [UT]{t.iso}', fontsize=20)
        plt.ylabel('Altitude [deg]', fontsize=20)
        plt.grid('both', linestyle='--', color='grey', alpha=0.5)
        plt.tight_layout()
        plt.axhline(y=altlimit*u.deg, linestyle='--', color='k', alpha=0.5)
        if not key_detect:
            plt.text(1, 75, 'Undetectable', fontsize = 15, c = 'r')
            plt.ylim(-40*u.deg, 90*u.deg)
            plt.savefig(f'{conf_param["savepath"]}undet_{ID}.png')
        else:
            plt.scatter((fir_time - t).jd*24, fir_alt, marker = '*', c= 'r')
            plt.text((fir_time - t).jd*24+0.3, fir_alt, 'First obstime = [UT] %.2d:%.2d:%.2d'%(fir_time.to_datetime().hour,fir_time.to_datetime().minute,fir_time.to_datetime().second), fontsize = 15, c = 'r')
            plt.ylim(0*u.deg, 90*u.deg)
            if key_detect == 2:
                plt.savefig(f'{conf_param["savepath"]}tm_det_{ID}.png')
            else:
                plt.savefig(f'{conf_param["savepath"]}det_{ID}.png')
    return fir_time.isot, key_detect
    
def handler(payload, root):
    
    #### Set required parameters from xml file ####
    pos2d = root.find('.//{*}Position2D')
    when = root.find(".//{*}Time")
    ID = str(root.find(".//Param[@name='TrigID']").attrib['value'])
    ra = float(pos2d.find('.//{*}C1').text)
    dec = float(pos2d.find('.//{*}C2').text)
    isotime = when.find(".//{*}ISOTime").text
    ###############################################
    
    fir_time, key_detect= check_detectability(ra,
                                              dec,
                                              conf_param['observatory'],
                                              isotime, check = True,
                                              frame = 'fk5',
                                              sunlimit = conf_param['sunlimit'],
                                              altlimit = conf_param['altlimit'],
                                              ID = ID)
    
    fir_time = Time(fir_time, format='isot',scale='utc').datetime
    
    if key_detect:
        start = '%d/%.2d/%.2d'%(fir_time.year,fir_time.month,fir_time.day)
        rtsmaker(targetcat=conf_param['catpath'],
                 observatory=conf_param['observatory'],
                 save_path=conf_param['savepath'],
                 start=start,
                 end=None,
                 headname='IMSNG',
                 altlimit=conf_param['altlimit'],
                 moonseperation=30.,
                 moonphase = True,
                 sunlimit=conf_param['sunlimit'],
                 numlimit=500)

        

#%%
    
    
    
    
    
    
    
    
    
    
    
#%%
import glob
path_test_alerts = sorted(glob.glob('/home/hhchoi1022/Desktop/Gitrepo/observatory/GCN_test/*/*.xml'))

#for xmlfile in path_test_alerts:
    xmlfile = path_test_alerts[0]
    try:
        payload = open(xmlfile, 'rb').read()
        root = lxml.etree.fromstring(payload)
        handler(payload, root)
    except:
        pass
    
