# Offered by Gu, Lim on Dec 2021  
# Modified by Hyeonho, Choi, 2022

#%%
def rtsmaker(observatory        = 'CBNUO',
             save_path          = './',
             obspath            = "/home/hhchoi1022/Desktop/Gitrepo/config/obs_location.txt",
             catpath            = './gecko_target.dat',
             start              = '2023/01/01',
             end                = '2023/01/03',
             headname           = 'IMSNG',
             altlimit           = 30.,
             moonseperation     = 40.,
             moonphase          = False,
             sunlimit           = -18,
             numlimit           = 500):
    import os, sys
    import pytz
    import jdcal
    import ephem
    import string
    import datetime
    import numpy as np
    import astropy.units as u
    from astropy.io import ascii
    import mskpy.observing as obs
    from astropy import units as u
    import astropy.coordinates as coord
    from astropy.coordinates import SkyCoord
    #------------------------------------------------------------#
    #    OBSERVATORY INFO.
    #------------------------------------------------------------#
    obsinfo          = ascii.read(obspath)
    obsname          = np.copy(obsinfo['name'])
    obsindex         = np.where(obsname == observatory)[0]
    obslat           = (np.copy(obsinfo['latitude(N+)'])[obsindex])[0]
    obslon           = (np.copy(obsinfo['longitude(E+)'])[obsindex])[0]
    obsalt           = (np.copy(obsinfo['altitude'])[obsindex])[0]
    obstz            = (np.copy(obsinfo['timezone'])[obsindex])[0]
    tz               = pytz.timezone(obstz)
    #------------------------------------------------------------#
    observ           = ephem.Observer()
    observ.lat       = str(obslat)
    observ.lon       = str(obslon)
    observ.elevation = obsalt
    observ.horizon   = str(sunlimit)
    #------------------------------------------------------------#
    #objects from catalog file
    tdata    = ascii.read(catpath)
    tdata.sort('ra')
    #pidx     = np.where( tdata['priority'] < 2 )[0]
    #tdata    = tdata[pidx]
    objname  = tdata['obj']
    ra       = tdata['ra']
    dec      = tdata['dec']
    prior    = tdata['priority']         
    dist     = tdata['dist']
    amaj     = tdata['maxaxis']
    bmin     = tdata['minaxis']
    RA       = coord.Angle(ra, unit = u.deg)
    DEC      = coord.Angle(dec, unit = u.deg)
    rad      = RA.degree
    radd     = RA.value
    decd     = DEC.degree
    decdd    = DEC.value

    #------------------------------------------------------------#
    #angular distance calculation
    def angsep(ra1deg, dec1deg, ra2deg, dec2deg) : 
        ra1rad   = ra1deg*(np.pi/180)
        dec1rad  = dec1deg*(np.pi/180)
        ra2rad   = ra2deg*(np.pi/180)
        dec2rad  = dec2deg*(np.pi/180)
        cos_a    = np.sin(dec1rad)*np.sin(dec2rad)+(np.cos(dec1rad)*np.cos(dec2rad)*np.cos(ra1rad-ra2rad))
        anglesep = np.arccos(cos_a)*180/np.pi
        return anglesep
    
    #------------------------------------------------------------#
    #dates to calculate
    fmt      = '%Y/%m/%d'
    if (start == None) & (end == None):
        startdt = datetime.datetime.today()
        enddt = startdt + datetime.timedelta(days = 1)
    if (start != None) & (end == None):
        startdt  = datetime.datetime.strptime(start, fmt)
        enddt = startdt + datetime.timedelta(days = 1)
    else:
        startdt  = datetime.datetime.strptime(start, fmt)
        enddt    = datetime.datetime.strptime(end, fmt)
    startmjd = (jdcal.gcal2jd(startdt.year, startdt.month, startdt.day))[1]
    endmjd   = (jdcal.gcal2jd(enddt.year, enddt.month, enddt.day))[1]
    
    for i in range(int(endmjd-startmjd)):
        onedaymjd   = startmjd+i
        filenameday = jdcal.jd2gcal(2400000.5, onedaymjd) # for filename
        oneday      = jdcal.jd2gcal(2400000.5, onedaymjd+1) # for real obs.
        onedaydt    = datetime.datetime(oneday[0], oneday[1], oneday[2], tzinfo=tz)
        dst         = tz.dst(onedaydt, is_dst=True)
        dst         = dst.seconds/3600

        onedayutc   = onedaydt.astimezone(pytz.utc)
        observ.date = onedayutc

        # Moon distance and information
        mcoord       = ephem.Moon()
        mcoord.compute(observ)
        minfo        = 'Moon ra, dec : '+str(mcoord.ra)+' '+str(mcoord.dec)+'\n'
        mphase       = ephem.Moon(observ.date)
        mphasestr    = 'Moon phase   : '+ "%.2f" % mphase.moon_phase +'\n'
        msep         = angsep(rad, decd, np.degrees(mcoord.ra), np.degrees(mcoord.dec))
        if moonphase:
            moonseperation = int(moonseperation/(1+np.power(1-mphase.moon_phase,2))) # Moonphase consideration when moonphase == True by Hyeonho Choi 
        
        #    SUNSET CALC.
        sunset      = observ.previous_setting(ephem.Sun())
        sunsettu    = ephem.Date.tuple(sunset)
        sunsetdt    = datetime.datetime(sunsettu[0],sunsettu[1],sunsettu[2],sunsettu[3],int(sunsettu[4]),tzinfo=pytz.utc)
        sunsetlocal = sunsetdt.astimezone(tz)
        sunsetstr   = str(sunlimit)+' deg sunset : '+str(sunsetlocal.hour)+':'+str(sunsetlocal.minute)+'\n'
        sunsethour  = sunsetlocal.hour+sunsetlocal.minute/60.+sunsetlocal.second/3600.
        #    SUNRISE CALC.
        sunrise      = observ.next_rising(ephem.Sun())
        sunrisetu    = ephem.Date.tuple(sunrise)
        sunrisedt    = datetime.datetime(sunrisetu[0],sunrisetu[1],sunrisetu[2],sunrisetu[3],int(sunrisetu[4]),tzinfo=pytz.utc)
        sunriselocal = sunrisedt.astimezone(tz)
        sunrisestr   = str(sunlimit)+' deg sunrise : '+str(sunriselocal.hour)+':'+str(sunriselocal.minute)+'\n'
        sunrisehour  = sunriselocal.hour+sunriselocal.minute/60.+sunriselocal.second/3600.

        stryear        = str(filenameday[0])
        strmonth       = str(filenameday[1])
        strday         = str(filenameday[2])
        if int(strmonth)    < 10 : strmonth = '0' + strmonth
        if int(strday)      < 10 : strday   = '0' + strday
        
        f = open(save_path+'/'+headname+'-'+stryear+strmonth+strday+"-rts_vis-"+observatory+".txt",'w')
        f.write('#\t'+str(observ.date)+' UTC & Day Time Saving +'+str(dst)+'\n')
        f.write('#\tObservatory\t= '+observatory+'\n')
        f.write('#\t'+sunsetstr)
        f.write('#\t'+sunrisestr)
        f.write('#\t'+minfo)
        f.write('#\t'+mphasestr)
        f.write('#\tMoon seperation = '+str(moonseperation)+'\n')
        f.write('#\tAltitude limit = '+str(altlimit)+'\n')
        f.write('#------------------------------------------------------- \n')
        f.write("name ra dec rise(LT) transit(LT) set(LT) moon_dist(deg) distance(Mpc) amaj['] bmin['] priority\n")
        numcount    = 0
        for n in range(len(rad)):
            param_rts    = dict(ra=rad[n],
                                dec=decd[n],
                                date=onedaydt,
                                lon=obslon,
                                lat=obslat,
                                tz=obstz,
                                limit=altlimit,
                                precision=1440)
            rtscal        = obs.rts(**param_rts)
            rt            = rtscal[0]
            tt            = rtscal[1]
            st            = rtscal[2]
            if rtscal[0]== None:
                pass
            elif sunrisehour < rtscal[0] < sunsethour and sunrisehour < rtscal[2] < sunsethour and sunrisehour < rtscal[1] < sunsethour: 
                pass
            elif msep[n] < moonseperation or msep[n] > 360-moonseperation:
                pass
            else:
                if numcount < numlimit:

                    nra = ra[n]
                    ndec = dec[n]
                    
                    rtp    = "%.2d" % int(rt)+':'+"%.2d" % int((rt-int(rt))*60)
                    ttp    = "%.2d" % int(tt)+':'+"%.2d" % int((tt-int(tt))*60)
                    stp    = "%.2d" % int(st)+':'+"%.2d" % int((st-int(st))*60)
                    vis = str(objname[n]) + '    ' + str(nra) + '    ' + str(ndec) + '    ' + rtp + '    ' + ttp + '    ' + stp + '    ' +str(int(msep[n])) + '     ' + str(int(dist[n]))+ '    '+str(amaj[n])+'   '+str(bmin[n])+'   '+str(prior[n])+'   '+'\n'
                    f.write(vis)
                    numcount+= 1
                else:
                    pass
        f.close()
#%%