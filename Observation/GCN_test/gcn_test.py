#!/usr/bin/env python
import pytz
import glob
import gcn
import lxml.etree
from astropy.time import Time
import numpy as np
import datetime
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.visualization import astropy_mpl_style, quantity_support
import matplotlib.pyplot as plt
plt.style.use(astropy_mpl_style)
quantity_support()
# Define your custom handler here.
# https://gcn.gsfc.nasa.gov/admin/gcn_swift3.txt
# https://github.com/lpsinger/pygcn
# https://gcn.gsfc.nasa.gov/vo_examples.html
@gcn.include_notice_types(
    gcn.notice_types.FERMI_GBM_FLT_POS,  # Fermi GBM localization (flight)
    gcn.notice_types.FERMI_GBM_GND_POS,  # Fermi GBM localization (ground)
    gcn.notice_types.FERMI_GBM_FIN_POS)  # Fermi GBM localization (final)

def handler(payload, root):

    trig = str(root.find(".//Param[@name='TrigID']").attrib['value'])
    #   Position
    pos2d = root.find('.//{*}Position2D')
    ra = float(pos2d.find('.//{*}C1').text)
    dec = float(pos2d.find('.//{*}C2').text)
    #	Target
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    radius = float(pos2d.find('.//{*}Error2Radius').text)
    gal_lat=root.find(".//Param[@name='Galactic_Lat']").attrib['value']

    #   Time
    when = root.find(".//{*}Time")
    dateobs = when.find(".//{*}ISOTime").text
    t = Time(dateobs, format='isot', scale='utc')

    #   Attribute
    not_GRB_label=root.find(".//Param[@name='Def_NOT_a_GRB']").attrib['value']
    test_lab=root.find(".//Param[@name='Test_Submission']").attrib['value']

    #   Catalog
    flt_cat_label=root.find(".//Param[@name='Target_in_Flt_Catalog']").attrib['value']
    gnd_cat_label=root.find(".//Param[@name='Target_in_Gnd_Catalog']").attrib['value']
    blk_cat_label=root.find(".//Param[@name='Target_in_Blk_Catalog']").attrib['value']

    #   Observatory
    timezone = 'Chile/Continental'
    tz = pytz.timezone(timezone)
    lat, lon, height =  -31.273333, -70.864444, 1122
    loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)

    t_local = pytz.utc.localize(t.datetime).astimezone(tz)  #   UTC --> Local time

    utcoffset = (t_local.utcoffset().total_seconds()*u.second).to(u.hour)
    midnight = t - utcoffset	#	local time
    # midnight = t
    
    #   Now + [0-24] hrs
    delta_midnight = np.linspace(0, 24, 200)*u.hour
    times = midnight + delta_midnight
    frame = AltAz(obstime=times, location=loc)
    c_altaz = c.transform_to(frame)

    #	Query Sun
    from astropy.coordinates import get_sun
    c_sun_altaz = get_sun(times).transform_to(frame)
    
    indx_sun = np.where(
        (c_sun_altaz.alt < -0*u.deg) &
        (c_altaz.alt > 30*u.deg)
    )
    ma = np.max(c_altaz.alt[indx_sun])    #   Meridian altitude

    #	Query Moon
    from astropy.coordinates import get_moon
    c_moon_altaz = get_moon(times).transform_to(frame)
    indx_moon, sep_moon, _ = c_altaz.match_to_catalog_sky(c_moon_altaz) #  Moon distance
    sep_moon_sun = sep_moon[indx_sun][sep_moon[indx_sun]>40*u.deg]  #   Moon distance (@night)

    if (ma > 30*u.deg) & (len(sep_moon_sun) > 0):
        #	Plot
        plt.close('all')
        plt.rcParams["figure.figsize"] = (15,6)
        plt.plot(delta_midnight, c_sun_altaz.alt, color='tomato', label='Sun', alpha=0.5)
        plt.plot(delta_midnight, c_moon_altaz.alt, color='dodgerblue', ls='--', label='Moon', alpha=0.5)
        plt.scatter(delta_midnight, c_altaz.alt, c=c_altaz.az, label=f'Swift_BAT_{trig}', lw=0, s=8, cmap='viridis')
        plt.fill_between(delta_midnight, 0, 90, c_sun_altaz.alt < -0)
        plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, c_sun_altaz.alt < -18*u.deg, color='0.65', zorder=0)

        plt.axhline(y=30*u.deg, linestyle='--', color='k', alpha=0.5)

        #	Setting
        plt.colorbar().set_label('Azimuth [deg]')
        plt.legend(loc='upper center', fontsize=16, framealpha=1.0)
        plt.xlim([np.min(delta_midnight), np.max(delta_midnight)])
        # plt.xticks((np.arange(13)*2-12)*u.hour)
        # plt.xticks(np.arange(0, 24+1)*u.hour)
        plt.xticks(np.arange(0, 24+1))
        plt.ylim(0*u.deg, 90*u.deg)
        # plt.xlabel('Hours from EDT Midnight', fontsize=20)
        plt.xlabel(f'Hours after {dateobs}', fontsize=20)
        plt.ylabel('Altitude [deg]', fontsize=20)
        plt.grid('both', linestyle='--', color='grey', alpha=0.5)
        plt.tight_layout()
        plt.tight_layout()
#%%
# Listen for VOEvents until killed with Control-C.
gcn.listen(handler=handler)
# path_test_alerts = sorted(glob.glob('/data1/GECKO/gcn/public_voe_examples/swift_test/*.xml'))
# for path_test_alert in path_test_alerts:
# path_test_alert = path_test_alerts[3]
# payload = open(path_test_alert, 'rb').read()
# root = lxml.etree.fromstring(payload)
handler(payload, root)
