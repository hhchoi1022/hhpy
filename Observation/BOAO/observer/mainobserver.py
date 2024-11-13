

#%%

from astropy.coordinates import EarthLocation, AltAz, SkyCoord, get_sun, get_moon
import astropy.units as u
from datetime import datetime
from astropy.time import Time
import pytz
from astroplan import Observer
from astropy.io import ascii
#%%
class mainObserver:
    '''
    This module define the basic information of the observatory.
    Most functions are originated from "Astroplan"
    See more info : https://astroplan.readthedocs.io
    '''
    def __init__(self, 
                 OBSERVER_LATITUDE : str = 36.1648, 
                 OBSERVER_LONGITUDE : str = 128.9766, 
                 OBSERVER_ELEVATION : str = 1162, 
                 OBSERVER_TIMEZONE : str = "Asia/Seoul",
                 OBSERVER_NAME : str = "Yunyi Choi & Hyeonho Choi",
                 OBSERVER_OBSERVATORY : str = "BOAO",
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
    
    ############ Time ############
    def localtime(self, 
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
    
    def siderialtime(self,
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
        
    def is_night(self,
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
    
    def tonight(self,
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
    
    def to_altaz(self,
                 ra : float or str or list,
                 dec : float or str or list,
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
        coord = self._observer.altaz(time, target = self.get_skycoord(ra, dec))
        return coord.altaz
    
    def to_radec(self,
                 alt : float,
                 az : float,
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
        coord = SkyCoord(alt = alt, az = az, frame ='altaz', location = self._earthlocation, obstime = time, unit = 'deg')
        return coord.icrs
    
    def risetime(self,
                 ra : float or str or list,
                 dec : float or str or list,
                 time : datetime or Time = None,
                 mode : str = 'next',
                 horizon : float = -18,
                 localtime : bool = True):
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
        if not isinstance(self.get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified')
        rise_time_utc = self._observer.target_rise_time(time, self.get_skycoord(ra, dec), which = mode, horizon = horizon*u.deg)
        rise_time_ltc = self.localtime(rise_time_utc)
        if localtime:
            return rise_time_ltc
        else:
            return rise_time_utc
    
    def settime(self,
                ra : float or str or list,
                dec : float or str or list,
                time : datetime or Time = None,
                mode : str = 'nearest',
                horizon : float = -18,
                localtime : bool = True):
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
        if not isinstance(self.get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified')
        rise_time_utc = self._observer.target_set_time(time, self.get_skycoord(ra, dec), which = mode, horizon = horizon*u.deg)
        rise_time_ltc = self.localtime(rise_time_utc)
        if localtime:
            return rise_time_ltc
        else:
            return rise_time_utc
    
    def meridiantime(self,
                     ra : float or str or list,
                     dec : float or str or list,         
                     time : datetime or Time = None,
                     mode : str = 'nearest',
                     localtime : bool = True):
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
        if not isinstance(self.get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified')
        rise_time_utc = self._observer.target_meridian_transit_time(time, self.get_skycoord(ra, dec), which = mode)
        rise_time_ltc = self.localtime(rise_time_utc)
        if localtime:
            return rise_time_ltc
        else:
            return rise_time_utc    
    
    def hourangle(self,
                  ra : float or str or list,
                  dec : float or str or list,   
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
        if not isinstance(self.get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified for hourangle')
        return self._observer.target_hour_angle(time, self.get_skycoord(ra, dec))                
    
    def parallactic_angle(self,
                          ra : float or str or list,
                          dec : float or str or list,   
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
        if not isinstance(self.get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified for hourangle')
        paralangle = self._observer.parallactic_angle(time = time, target = self.get_skycoord(ra, dec)).deg
        #if paralangle < 0:
        #    paralangle += 360
        print('Parallactic angle = %.1f'%paralangle)
        return paralangle

    ############ Sun ############
    def sun_radec(self,
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
    
    def sun_altaz(self,
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
    
    def sun_risetime(self,
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
    
    def sun_settime(self,
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
    def moon_radec(self,
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
    
    def moon_altaz(self,
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
            
        return self.to_altaz(radec = get_moon(time), time = time)

    def moon_risetime(self,
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
    
    def moon_settime(self,
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
    
    def moon_phase(self,
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
    
    def get_reference_star(self,
                           ra : float or str,
                           dec : float or str,
                           time : datetime or Time = None,
                           idx_from_smallest : int = 0,
                           refstar_path : str = './StdStars.txt',
                           refstar_ra_key : str = 'ra',
                           refstar_dec_key : str = 'dec',
                           sort_by : str = 'separation',
                           altitude_cut : float = 30,
                           consider_meridian_flip: bool = True,
                           meridian_flip_angle_tolerance: float = 15):
        if time == None:
            time = Time.now()
        if not isinstance(time, Time):
            time = Time(time)
            
        target_coord = self.get_skycoord(ra, dec)
        target_altaz = self.to_altaz(target_coord.ra.deg, target_coord.dec.deg, time = time)
        refstar_tbl = ascii.read(refstar_path)
        refstar_coord = self.get_skycoord(refstar_tbl[refstar_ra_key], refstar_tbl[refstar_dec_key])
        refstar_altaz = self.to_altaz(refstar_coord.ra.deg, refstar_coord.dec.deg, time = time)
        separation = target_coord.separation(refstar_coord).deg
        refstar_tbl['alt'] = refstar_altaz.alt.deg
        refstar_tbl['az'] = refstar_altaz.az.deg
        refstar_tbl['alt_diff'] = abs(target_altaz.alt.deg - refstar_altaz.alt.deg)
        refstar_tbl['az_diff'] = abs(target_altaz.az.deg - refstar_altaz.az.deg)
        refstar_tbl['separation'] = separation
        
        refstar_tbl = refstar_tbl[refstar_tbl['alt'] > altitude_cut]
        if len(refstar_tbl) == 0:
            raise ValueError('No reference star is found due to altitude criteria')
        
        if sort_by.upper() == 'SEPARATION':
            refstar_tbl.sort(['separation'])
        elif sort_by.upper() == 'ALTITUDE':
            refstar_tbl.sort('alt_diff')
        elif sort_by.upper() == 'AZIMUTH':
            refstar_tbl.sort('az_diff')
        
        if consider_meridian_flip:
            direction_target = target_altaz.az.deg > 180
            if direction_target:
                filter_refstar  = refstar_altaz.az.deg > 180
            else:
                filter_refstar  = refstar_altaz.az.deg < 180 - meridian_flip_angle_tolerance
            refstar_tbl = refstar_tbl[filter_refstar]
            if len(refstar_tbl) == 0:
                raise ValueError('No reference star is found due to meridian flip criteria')
            
        refstar = refstar_tbl[idx_from_smallest]
        print('Target Alt = %.1f, Az = %.1f'%(target_altaz.alt.deg, target_altaz.az.deg))
        print('Ref Alt = %.1f, Az = %.1f'%(refstar['alt'], refstar['az']))
        print('Ref Separation = %.1f'%refstar['separation'])
        return refstar
        
    def get_skycoord(self,
                     ra : str or float or list,
                     dec: str or float or list):
        """
        Parameters
        ==========
        ra : str | float | list = Right Ascension, if str, it should be in hms format (e.g., "10:20:30"), if float, it should be in decimal degrees
        dec : str | float | list = Declination, if str, it should be in dms format (e.g., "+20:30:40"), if float, it should be in decimal degrees
        
        Return
        ======
        coord : SkyCoord = SkyCoord object
        """
        
        # Check if RA and Dec are given as strings (like "10:20:30")
        if isinstance(ra, str) and isinstance(dec, str):
            # Interpret as sexagesimal format (e.g., "10:20:30", "+20:30:40")
            coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        elif isinstance(ra, (float, int)) and isinstance(dec, (float, int)):
            # Interpret as decimal degrees
            coord = SkyCoord(ra * u.deg, dec * u.deg)
        else:
            # Interpret as list of RA and Dec
            if len(ra) != len(dec):
                raise ValueError("RA and Dec should have the same length")
            if isinstance(ra[0], str) and isinstance(dec[0], str):
                coord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg))
            elif isinstance(ra[0], (float, int)) and isinstance(dec[0], (float, int)):
                coord = SkyCoord(ra = ra*u.deg, dec = dec*u.deg)
            else:
                raise ValueError("Unsupported RA and Dec format")
        return coord
    
    def get_object_radec(self,
                         objname: str,
                         obj_path : str = './ObservableTargets1.txt',
                         obj_key : str = 'Name',
                         ra_key : str = 'ra',
                         dec_key : str = 'dec'
                         ):
        """
        Parameters
        ==========
        objname : str = Name of the object
        obj_path : str = Path to the object list
        
        Return
        ======
        ra : float = RA of the object
        dec : float = DEC of the object
        """
        obj_tbl = ascii.read(obj_path)
        if not objname in obj_tbl[obj_key]:
            raise ValueError('No object is found \n Registered object names are \n %s'%obj_tbl[obj_key])
        else:
            target = obj_tbl[obj_tbl[obj_key] == objname]      
            return target[ra_key][0], target[dec_key][0]
        
# %% fractical use
if __name__ == '__main__':
    observer = mainObserver()
    ra, dec = observer.get_object_radec('J2139+2929')
    refstar = observer.get_reference_star(ra,
                                          dec,
                                          time= None,
                                          idx_from_smallest = 0,
                                          refstar_path = './StdStars.txt',
                                          refstar_ra_key = 'ra',
                                          refstar_dec_key = 'dec',
                                          sort_by = 'separation',
                                          altitude_cut = -60,
                                          consider_meridian_flip = True,
                                          meridian_flip_angle_tolerance = 15)
    parallel_angle = observer.parallactic_angle(ra, dec)
#%%
    
