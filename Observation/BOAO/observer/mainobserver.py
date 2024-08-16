

#%%

from astropy.coordinates import EarthLocation, AltAz, SkyCoord, get_sun, get_moon
import astropy.units as u
from datetime import datetime
from astropy.time import Time
import pytz
from astroplan import Observer

#%%
class mainObserver:
    '''
    This module define the basic information of the observatory.
    Most functions are originated from "Astroplan"
    See more info : https://astroplan.readthedocs.io
    '''
    def __init__(self, 
                 OBSERVER_LATITUDE : str, 
                 OBSERVER_LONGITUDE : str, 
                 OBSERVER_ELEVATION : str, 
                 OBSERVER_TIMEZONE : str,
                 OBSERVER_NAME : str = None,
                 OBSERVER_OBSERVATORY : str = None,
                 **kwargs
                 ):
        """
        Parameters
        ==========
        OBSERVER_LATITUDE : str = Earth latitude
        OBSERVER_LONGITUDE : str = Earth East longitude
        OBSERVER_ELEVATION : str = Height above the Earth surface (reference ellipsoid = 'WGC84)
        OBSERVER_TIMEZONE : str = Timezone information, for all list of timezone, <> pytz.all_timezones
        OBSERVER_NAME : str = Name of the observer(None)
        OBSERVER_OBSERVATORY : str = Name of the observatory(None)
        
        Functions
        ======
        
        """
        self._latitude = None
        self._longitude = None
        self._elevation = 0
        self._timezone = None
        self._name = None
        self._observatory = None
        
        if OBSERVER_LATITUDE is not None:
            self._latitude = float(OBSERVER_LATITUDE)*u.deg
        if OBSERVER_LONGITUDE is not None:
            self._longitude = float(OBSERVER_LONGITUDE)*u.deg
        if OBSERVER_ELEVATION is not None:        
            self._elevation = float(OBSERVER_ELEVATION)*u.m
        if OBSERVER_NAME is not None:        
            self._name = OBSERVER_NAME
        if OBSERVER_OBSERVATORY is not None:        
            self._observatory= OBSERVER_OBSERVATORY
        if OBSERVER_TIMEZONE is not None:
            self._timezone = pytz.timezone(OBSERVER_TIMEZONE)
        if (OBSERVER_LATITUDE is not None) & (OBSERVER_LONGITUDE is not None):
            self._earthlocation = EarthLocation.from_geodetic(lat=self._latitude, lon=self._longitude, height=self._elevation)
            self._observer = Observer(location = self._earthlocation, name = self._observatory, timezone = self._timezone)
    ############ Site info ############
    def obs_latitude(self):
        return self._latitude
    def obs_longitude(self):
        return self._longitude
    def obs_elevation(self):
        return self._elevation
    def obs_observername(self):
        return self._name
    def obs_observatoryname(self):
        return self._observatory
    def obs_observer(self):
        return self._observer
    
    ############ Time ############
    def obs_localtime(self, 
                      utctime : datetime = None):
        """
        Parameters
        ==========
        utctime : datetime or Time = Time(default = Now)
        
        Return
        ======
        localtime : datetime = local time for the given site
        """
        if utctime == None:
            utctime = datetime.utcnow()
        localtime = pytz.utc.localize(utctime).astimezone(self._timezone)
        return localtime
    
    def obs_siderialtime(self,
                         time : datetime or Time = None,
                         mode : str = 'mean'): # Or apparent 
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        mode : str = mean or apparent, mode to calculate sidereal time(default = mean)
        
        Return
        ======
        siderialtime : Longitude = siderialtime time for the given site
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.local_sidereal_time(time, kind = mode)
    
    def obs_now(self):
        return Time.now()
    
    def obs_is_night(self,
                     time : datetime or Time = None):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        is_night : bool = True if the time is night
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.is_night(time, horizon = -18*u.deg)
    
    def obs_tonight(self,
                    time : datetime or Time = None,
                    horizon = -18):
        """
        Parameters
        ==========
        time : datetime or Time = Time to calculate tonight(default = Now)
        horizon : float = Degrees above/below actual horizon to use for calculating rise/set times
        
        Return
        ======
        tonight : tuple = tuple[0] for night start time, tuple[1] for night end time
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.tonight(time, horizon = horizon*u.deg)

    ############ Target ############
    
    def obs_to_altaz(self,
                     radec : SkyCoord,
                     time : datetime or Time = None):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to convert AltAz
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        coordinate : SkyCoord = AltAz coordinate of the target
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.altaz(time, target = radec)
    
    def obs_to_radec(self,
                     altaz : SkyCoord,
                     time : datetime or Time = None):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(AltAz) coordinates to convert RADec
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        coordinate : SkyCoord = RADec coordinate of the target
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        coord = SkyCoord(alt = altaz.alt, az = altaz.az, frame ='altaz', location = self._earthlocation, obstime = time)
        return coord.icrs
    
    def obs_risetime(self,
                     radec : SkyCoord,
                     time : datetime or Time = None,
                     mode : str = 'next',
                     horizon : float = -18):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate hourangle
        time : datetime or Time = Time(default = Now)
        mode : {‘next’, ‘previous’, ‘nearest’}(default = next)
        horizon : float = horizon angle in degree(default = -18)
        
        Return
        ======
        risetime : Time = Target risetime
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        if not isinstance(radec, SkyCoord):
            raise ValueError('No target is specified')
        return self._observer.target_rise_time(time, radec, which = mode, horizon = horizon*u.deg)
    
    def obs_settime(self,
                    radec : SkyCoord,
                    time : datetime or Time = None,
                    mode : str = 'nearest',
                    horizon : float = -18):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate hourangle
        time : datetime or Time = Time(default = Now)
        mode : {‘next’, ‘previous’, ‘nearest’}(default = nearest)
        horizon : float = horizon angle in degree(default = -18)
        
        Return
        ======
        settime : Time = Target settime
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        if not isinstance(radec, SkyCoord):
            raise ValueError('No target is specified')
        return self._observer.target_set_time(time, radec, which = mode, horizon = horizon*u.deg)
    
    def obs_meridiantime(self,
                         radec : SkyCoord,
                         time : datetime or Time = None,
                         mode : str = 'nearest'):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate hourangle
        time : datetime or Time = Time(default = Now)
        mode : {‘next’, ‘previous’, ‘nearest’}(default = nearest)
        
        Return
        ======
        meridiantime : Time = Target meridian transit time
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        if not isinstance(radec, SkyCoord):
            raise ValueError('No target is specified')
        return self._observer.target_meridian_transit_time(time, radec, which = mode)
    
    def obs_hourangle(self,
                      radec : SkyCoord,
                      time : datetime or Time = None):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate hourangle
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        hourangle : Longitude = Hourangle of the target, + For west(pass the meridian)
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        if not isinstance(radec, SkyCoord):
            raise ValueError('No target is specified for hourangle')
        return self._observer.target_hour_angle(time, radec)                
    
    def obs_parallactic_angle(self,
                              radec : SkyCoord,
                              time : datetime or Time = None):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate parallactic angle
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        parallangle : Angle = parallactic_angle of the target
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        if not isinstance(radec, SkyCoord):
            raise ValueError('No target is specified for hourangle')
        paralangle = self._observer.parallactic_angle(time = time, target = radec).deg
        #if paralangle < 0:
        #    paralangle += 360
        return paralangle

    ############ Sun ############
    def obs_sun_radec(self,
                      time : datetime or Time = None):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        sunradec : SkyCoord = RADec coordinate of the Sun
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return get_sun(time)
    
    def obs_sun_altaz(self,
                      time : datetime or Time = None):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        sunaltaz : SkyCoord = AltAz coordinate of the Sun
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.sun_altaz(time)
    
    def obs_sun_risetime(self,
                         time : datetime or Time = None,
                         mode = 'nearest',
                         horizon = -18):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        mode : {‘next’, ‘previous’, ‘nearest’}
        
        
        Return
        ======
        sunrisetime : Time = Sunrise time 
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.sun_rise_time(time, which = mode, horizon = horizon * u.deg)
    
    def obs_sun_settime(self,
                        time : datetime or Time = None,
                        mode = 'nearest',
                        horizon = -18):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        mode : str = mode {‘next’, ‘previous’, ‘nearest’}
        
        
        Return
        ======
        sunsettime : Time = Sunset time 
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.sun_set_time(time, which = mode, horizon = horizon * u.deg)
    
    ############ Moon ############
    def obs_moon_radec(self,
                       time : datetime or Time = None):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        moonradec : SkyCoord = RADec coordinate of the Moon
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return get_moon(time)
    
    def obs_moon_altaz(self,
                       time : datetime or Time = None):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        moonaltaz : SkyCoord = AltAz coordinate of the Moon
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
            
        return self.obs_to_altaz(radec = get_moon(time), time = time)

    def obs_moon_risetime(self,
                          time : datetime or Time = None,
                          mode = 'nearest',
                          horizon = -18):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        mode : str = mode {‘next’, ‘previous’, ‘nearest’}(default = nearest)
        horizon : float = horizon angle in degree(default = -18)
        
        Return
        ======
        risetime : Time = Moon rise time above the horizon
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.moon_rise_time(time, which = mode, horizon = horizon * u.deg)
    
    def obs_moon_settime(self,
                         time : datetime or Time = None,
                         mode = 'nearest',
                         horizon = -18):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        mode : str = mode {‘next’, ‘previous’, ‘nearest’}(default = nearest)
        horizon : float = horizon angle in degree(default = -18)
        
        
        Return
        ======
        settime : Time = Moon Set time above the horizon
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.moon_set_time(time, which = mode, horizon = horizon * u.deg)
    
    def obs_moon_phase(self,
                       time : datetime or Time = None):
        """
        Parameters
        ==========
        time : datetime or Time = Time(default = Now)
        
        Return
        ======
        phase : float = 0 is new, 1 is full
        """
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
        return self._observer.moon_illumination(time)
# %% fractical use
if __name__ == '__main__':
    ra = 31.5666
    dec = -0.2914
    import json
    with open('Observer_BOAO.config', 'r') as f:
        config = json.load(f)
    observer = mainObserver(**config)
    target = SkyCoord(ra, dec, unit = 'deg')
    paralangle = observer.obs_parallactic_angle(target, time = Time.now())
    hourangle = observer.obs_hourangle(target, time = Time.now()).hour
    print(f'Parallactic Angle = %.1f, Hour Angle = %.1f'%(paralangle, hourangle))
    t = Time(observer.obs_meridiantime(target), format = 'jd')
    maxaltaz = observer.obs_to_altaz(target, t)
    print('meridian time(UT) : %d/%d %d:%d:%d' %(t.datetime.month, t.datetime.day, t.datetime.hour, t.datetime.minute, t.datetime.second))
    print('Alt = %.1f, Az at meridian = %.1f'%(maxaltaz.az.value, maxaltaz.alt.value))
    #print(observer.obs_parallactic_angle(target, time = timelist))
#%%
