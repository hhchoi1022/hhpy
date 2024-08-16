#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 22:08:40 2022

@author: hhchoi1022
"""
#%%
conf_param = dict(
                  # Configuration path setting 
                  obspath          = "/home/hhchoi1022/Desktop/Gitrepo/config/obs_location.txt",
                  savepath         = '/home/hhchoi1022/ACPscript/ToO/', # path for saving the scripts  
                  # rts & ACP script setting 
                  catpath          = '/home/hhchoi1022/Desktop/Gitrepo/config/alltarget_KCT_220904.txt',
                  observatory      = 'KCT',
                  altlimit         = 35,
                  sunlimit         = -18,
                  moonseperation   = 30,
                  moonphase        = True,
                  overhead_ratio   = 1.1,
                  select_mode      = 'altitude',
                  dark             = True,
                  bias             = True,
                  flat             = True,
                  frame            = 'fk5',
                  #ToO observation setting 
                  ID               = 'GRB221009A',
                  importance       = 0.95,
                  ra               = 288.263,
                  dec              = +19.803,
                  exptime          = '120,120',
                  counts           = '10,10',
                  filter           = 'i,z',
                  binning          = '1,1',
                  )
#%%
import os
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
import datetime
from rtsmaker import rtsmaker
from ACPscript_v2 import  ACP_scriptmaker
from xml.etree.ElementTree import ElementTree
plt.style.use(astropy_mpl_style)
quantity_support()

############# ADD NOTICE SOURCES HERE #############
##list = https://gcn.gsfc.nasa.gov/filtering.html##

@gcn.include_notice_types(
    gcn.notice_types.FERMI_GBM_FLT_POS,  # Fermi GBM localization (flight)
    gcn.notice_types.FERMI_GBM_GND_POS,  # Fermi GBM localization (ground)
    gcn.notice_types.FERMI_GBM_FIN_POS,  # Fermi GBM localization (final)
    gcn.notice_types.LVC_PRELIMINARY, # LIGO, VIRGO
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE,
    #gcn.notice_types.LVC_RETRACTION,
    gcn.notice_types.SWIFT_BAT_GRB_POS_ACK
    #gcn.notice_types.SWIFT_XRT_POSITION,
    #gcn.notice_types.SWIFT_UVOT_POS
    )
############# ADD NOTICE SOURCES END ##############

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

def save_xml(root, filename):
    tree = ElementTree(root)
    tree.write(f'{conf_param["savepath"]}'+filename+'.xml')
    

def check_detectability(ra, dec, observatory, time, importance, graph = True, frame = 'fk5', sunlimit = -18, altlimit = 30, ID = None):
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
    ra, dec = to_skycoord(ra, dec).ra.value, to_skycoord(ra,dec).dec.value
    
    t = Time(time, format='isot', scale='utc')
    #t_now = datetime.datetime.utcnow()
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
    moondist_tar = np.median(sep_moon)
    mphase       = ephem.Moon(t.iso)
    # Object detectability
    if len(indx_obj[0])==0:
        fir_time = times[0]
        fir_alt = 0
        key_detect = 0
    
    else:
        fir_time = np.min(times[indx_obj])
        fir_alt = c_altaz.alt[indx_obj[0][np.argmin(times[indx_obj])]]
        if np.abs(t.datetime-fir_time.datetime).seconds/3600< 1:
            key_detect = 1
        else:
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
        plt.text(17, 85, f'Moon seperation = {int(moondist_tar.degree)}', fontsize = 15)
        plt.text(17, 80, 'Moon phase = %.2f' %mphase.moon_phase, fontsize = 15)
        plt.text(17, 75, 'Importance = %.2f' %importance, fontsize = 15)
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
            plt.text(17, 70, 'Undetectable', fontsize = 15, c = 'r')
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

def IDmaker():
    now = datetime.datetime.now()
    ID = '%.2d%.2d%.2d'%(now.year%2000,now.month,now.day)
    return ID

def main(ra, dec,
         isotime,
         ID,
         importance,
         exptime,
         counts,
         filter_,
         binning
         ):
    if importance > 0.9:
        ###############################################
        fir_time, key_detect= check_detectability(ra,
                                                  dec,
                                                  conf_param['observatory'],
                                                  isotime, 
                                                  graph = True,
                                                  frame = 'fk5',
                                                  sunlimit = conf_param['sunlimit'],
                                                  altlimit = conf_param['altlimit'],
                                                  ID = ID,
                                                  importance = importance)
        if type(ra) == float:
            coord = SkyCoord(ra, dec, unit = u.deg)
            ra = '%.2d:%.2d:%06.3f'%(coord.ra.hms.h,coord.ra.hms.m,coord.ra.hms.s)
            dec = '%.2d:%.2d:%05.2f'%(coord.dec.dms.d,np.abs(coord.dec.dms.m),np.abs(coord.dec.dms.s))
        fir_time = Time(fir_time, format='isot',scale='utc').datetime
        if fir_time.hour < 12: # When first detectable time pass the night
            fir_night = fir_time - datetime.timedelta(days=--1)
        else:
            fir_night = fir_time
        if key_detect == 1:
            start = '%d/%.2d/%.2d'%(fir_night.year,fir_night.month,fir_night.day)
            rtsfile = rtsmaker(
                      catpath=conf_param['catpath'],
                      observatory=conf_param['observatory'],
                      save_path=conf_param['savepath'],
                      obspath = conf_param['obspath'],
                      start=start,
                      end=None,
                      headname='IMSNG',
                      altlimit=conf_param['altlimit'],
                      moonseperation=conf_param['moonseperation'],
                      moonphase = True,
                      sunlimit=conf_param['sunlimit'],
                      numlimit=500
                      )
            scriptmaker = ACP_scriptmaker(rtspath = rtsfile,
                                          observatory=conf_param['observatory'],
                                          overhead_ratio = conf_param['overhead_ratio'],
                                          select_mode = conf_param['select_mode'],
                                          dark = conf_param['dark'],
                                          bias = conf_param['bias'],
                                          flat = conf_param['flat'])
            _ = scriptmaker.add_ToO(ID = ID, ra = ra, dec = dec, obstime = fir_time, exptime = exptime, counts = counts, filter_ = filter_, binning = binning)
            ToOfilename = os.path.dirname(scriptmaker.scriptpath)+'/ToO_'+os.path.basename(scriptmaker.scriptpath)
            scriptmaker.make_ACPscript(filename = ToOfilename)
            scriptmaker.savefigpath = ToOfilename[:-3]+'png'
            scriptmaker.show_graph()
            scriptmaker.cut_ToOscript(scriptpath = ToOfilename)    
        if key_detect == 2:
            start = '%d/%.2d/%.2d'%(fir_night.year,fir_night.month,fir_night.day)
            rtsfile = rtsmaker(
                      catpath=conf_param['catpath'],
                      observatory=conf_param['observatory'],
                      save_path=conf_param['savepath'],
                      obspath = conf_param['obspath'],
                      start=start,
                      end=None,
                      headname='IMSNG',
                      altlimit=conf_param['altlimit'],
                      moonseperation=conf_param['moonseperation'],
                      moonphase = True,
                      sunlimit=conf_param['sunlimit'],
                      numlimit=500
                      )
            scriptmaker = ACP_scriptmaker(rtspath = rtsfile,
                                          observatory=conf_param['observatory'],
                                          overhead_ratio = conf_param['overhead_ratio'],
                                          select_mode = conf_param['select_mode'],
                                          dark = conf_param['dark'],
                                          bias = conf_param['bias'],
                                          flat = conf_param['flat'])
            _ = scriptmaker.add_ToO(ID = ID, ra = ra, dec = dec, obstime = fir_time, exptime = exptime, counts = counts, filter_ = filter_, binning = binning)
            ToOfilename = os.path.dirname(scriptmaker.scriptpath)+'/ToO_'+os.path.basename(scriptmaker.scriptpath)
            scriptmaker.savefigpath = ToOfilename[:-3]+'png'
            scriptmaker.show_graph()
            scriptmaker.make_ACPscript(filename = os.path.dirname(scriptmaker.scriptpath)+'/ToO_tm_'+os.path.basename(scriptmaker.scriptpath))
        return key_detect

def handler(payload, root):
    
    ##### Set required parameters from xml file #####
    pos2d = root.find('.//{*}Position2D')
    when = root.find(".//{*}Time")
    ra = float(pos2d.find('.//{*}C1').text)
    dec = float(pos2d.find('.//{*}C2').text)
    isotime = when.find(".//{*}ISOTime").text
    time = datetime.datetime.strptime(isotime,'%Y-%m-%dT%H:%M:%S.%f')
    isotimeID = '{:02}{:02}{:02}'.format(time.year%2000, time.month, time.day)
    ID = 'ToO'+isotimeID+'_'+str(ra)+'_'+str(dec)
    why = root.find(".//{*}Why")
    importance = float(why.items()[0][1])
    key_detect = main(ra = ra, dec = dec, isotime = isotime, ID = ID, importance = importance, exptime = conf_param['exptime'], counts = conf_param['counts'], filter_ = conf_param['filter'], binning = conf_param['binning'])
    if (key_detect == 1) | (key_detect == 2):
        save_xml(root,filename = ID+'.xml')
#%%


if __name__ == '__main__':
    gcn.listen(handler=handler)
        
#%% Test
'''
import glob
path_test_alerts = sorted(glob.glob('/home/hhchoi1022/Desktop/Gitrepo/observatory/GCN_test/swift/*.xml'))

for xmlfile in path_test_alerts:
    try:
        payload = open(xmlfile, 'rb').read()
        root = lxml.etree.fromstring(payload)
        handler(payload, root)
    except:
        pass
'''
#%% Test for ToO

from HHsupport_phot import to_skycoord
isotime = datetime.datetime.utcnow().isoformat()
ra_deg, dec_deg = to_skycoord(conf_param['ra'], conf_param['dec']).ra.value, to_skycoord(conf_param['ra'], conf_param['dec']).dec.value
main(ra = ra_deg, dec = dec_deg, isotime = isotime, ID = conf_param['ID'], importance=  conf_param['importance'], exptime = conf_param['exptime'], counts = conf_param['counts'], filter_ = conf_param['filter'], binning = conf_param['binning'])


# %%
