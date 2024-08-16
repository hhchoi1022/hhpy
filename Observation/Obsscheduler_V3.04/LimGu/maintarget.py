#%%
# Other modules
from astroplan import FixedTarget, is_event_observable, is_observable, is_always_observable
from astroplan import AltitudeConstraint, AirmassConstraint, MoonSeparationConstraint, GalacticLatitudeConstraint, AtNightConstraint
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib import cm
# TCSpy modules
from mainobserver import mainObserver
from mainconfig import mainConfig
#%%
class mainTarget(mainConfig):
    """
    Parameters
    ----------
    1. observer : mainObserver
        An instance of mainObserver representing the observer.
    2. target_ra : float, optional
        The right ascension of the target, in hours.
    3. target_dec : float, optional
        The declination of the target, in degrees.
    4. target_alt : float, optional
        The altitude of the target, in degrees.
    5. target_az : float, optional
        The azimuth of the target, in degrees.
    6. target_name : str, optional
        The name of the target.

    Methods
    -------
    1. get_status() -> dict
        Returns a dictionary with information about the current status of the target.
    2. is_event_observable(utctimes: datetime or Time = None) -> bool
        Determines whether the target is observable at the specified time or at the current time.
    3. altaz(utctimes: datetime or Time = None) -> SkyCoord
        Calculate the alt-az coordinates of the target for a given time(s) in UTC.
    4. risetime(utctime: datetime or Time = None, mode: str = 'next', horizon: float = 30) -> Time
        Calculate the next rise time of the target as seen by the observer.
    5. settime(utctime: datetime or Time = None, mode: str = 'nearest', horizon: float = 30) -> Time
        Calculate the time when the target sets below the horizon.
    6. meridiantime(utctime: datetime or Time = None, mode: str = 'nearest') -> Time
        Calculate the time at which the target passes through the observer's meridian.
    7. hourangle(utctimes: datetime or Time = None) -> Angle
        Calculate the hour angle of the target(s) at the specified time(s).
    8. staralt(utctime : datetime or Time or np.array = None)
        Creates a plot of the altitude and azimuth of a celestial object.
    """
    
    def __init__(self,
                 name_telescope : str,
                 name_project : str,
                 observer : mainObserver,
                 target_ra : float = None,
                 target_dec : float = None,
                 target_alt : float = None,
                 target_az : float = None,
                 target_name : str = '',
                 **kwargs):
        
        super().__init__(name_telescope = name_telescope, name_project = name_project)
        self._observer = observer
        self._astroplan_observer = observer.info['observer']
        self._constraints = self._get_constraints(**self.config)
        self.ra = target_ra
        self.dec = target_dec
        self.alt = target_alt
        self.az = target_az
        
        if (isinstance(target_ra, type(None))) & (isinstance(target_dec, type(None))):
            self.alt = target_alt
            self.az = target_az
            self.coordinate = self._get_coordinate_altaz(alt = target_alt, az = target_az)
            self.target_astroplan = self._get_target(self.coordinate, target_name)
     
        else:
            self.coordinate = self._get_coordinate_radec(ra = target_ra, dec = target_dec)
            self.target_astroplan = self._get_target(self.coordinate, target_name)
            altaz = self.altaz()
            #self.ra = self.coordinate.ra.hour
            self.ra = self.coordinate.ra.deg
            self.dec = self.coordinate.dec.deg
            self.alt = altaz.alt.value
            self.az = altaz.az.value
            
        self.name = target_name
        #self.status = self.get_status()
        
    #def __repr__(self):
    #    return f'Target[name={self.name}, ra[deg]={round(self.ra,4)}, dec[deg]={round(self.dec,4)}]'

    def get_status(self):
        """
        Returns a dictionary with information about the current status of the target.
        
        Return
        ======
        targetinfo : dict
            A dictionary containing the following fields:
                - update_time: the current time in ISO format.
                - jd : the current time in JD format
                - name: the name of the target.
                - ra: the right ascension of the target.
                - dec: the declination of the target.
                - alt: the altitude of the target in degrees.
                - az: the azimuth of the target in degrees.
                - coordtype : the coordinate type defined ('radec' or 'altaz')
                - hourangle: the hour angle of the target in degrees.
                - is_event_observable: a boolean indicating whether the target is currently observable.
        """
        if self.coordinate.frame.name == 'altaz':
            targetinfo = dict()
            targetinfo['update_time'] = Time.now().isot
            targetinfo['jd'] = np.round(Time.now().jd,5)
            targetinfo['name'] = self.name
            targetinfo['ra'] = None
            targetinfo['dec'] = None
            targetinfo['alt'] = np.round(self.alt,3)
            targetinfo['az'] = np.round(self.az,3)
            targetinfo['coordtype'] = 'altaz'
            targetinfo['hourangle'] = None
            targetinfo['is_event_observable'] = None
        else:
            targetinfo = dict()
            targetinfo['update_time'] = Time.now().isot
            targetinfo['jd'] = np.round(Time.now().jd,5)
            targetinfo['name'] = self.name
            targetinfo['ra'] = np.round(self.ra,5)
            targetinfo['dec'] = np.round(self.dec,5)
            targetinfo['alt'] = np.round(self.alt,3)
            targetinfo['az'] = np.round(self.az,3)
            targetinfo['coordtype'] = 'radec'
            targetinfo['hourangle'] = np.round(self.hourangle().value,3)
            targetinfo['is_event_observable'] = self.is_event_observable()
        return targetinfo

    def is_ever_observable(self,
                           utctime : datetime or Time = None) -> bool:
        """
        Determines whether the target is observable during the specified time

        Parameters
        ----------
        1. utctimes : datetime or Time, optional
            The time at which to check observability. Defaults to the current time.
            
        Returns
        -------
        bool
            True if the target is observable, False otherwise.
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        tonight = self._observer.tonight(utctime)
        starttime = tonight[0]
        endtime = tonight[1]
        time_range = [starttime, endtime]
        return is_observable(constraints = self._constraints, observer = self._astroplan_observer, targets = self.target_astroplan, time_range = time_range)
    

    def is_always_observable(self,
                             utctimes : datetime or Time or np.array = None) -> bool:
        """
        Determines whether the target is always observable during the specified time

        Parameters
        ----------
        1. utctimes : datetime or Time, optional
            The time at which to check observability. Defaults to the current time.
            
        Returns
        -------
        bool
            True if the target is observable, False otherwise.
        """
        if utctimes is None:
            utctimes = Time.now()
        if not isinstance(utctimes, Time):
            utctimes = Time(utctimes)
        return is_always_observable(constraints = self._constraints, observer = self._astroplan_observer, targets = self.target_astroplan, times = utctimes)
    
    def is_event_observable(self,
                            utctimes : datetime or Time or np.array = None) -> bool:
        """
        Determines whether the target is observable at the specified time or at the current time.

        Parameters
        ----------
        1. utctimes : datetime or Time, optional
            The time at which to check observability. Defaults to the current time.
            
        Returns
        -------
        bool
            True if the target is observable, False otherwise.
        """
        if utctimes is None:
            utctimes = Time.now()
        if not isinstance(utctimes, Time):
            utctimes = Time(utctimes)
        return is_event_observable(constraints = self._constraints, observer = self._astroplan_observer, target = self.target_astroplan, times = utctimes)
    
    def altaz(self,
              utctimes : datetime or Time or np.array = None) -> SkyCoord:
        """
        Calculate the alt-az coordinates of the target for a given time(s) in UTC.

        Parameters
        ==========
        1. utctimes : datetime or Time, optional
            The time(s) to calculate the alt-az coordinates for, in UTC. If not provided, the current time will be used. 

        Returns
        =======
        1. SkyCoord
            The alt-az coordinates of the target at the specified time(s).
        """
        if utctimes is None:
            utctimes = Time.now()
        if not isinstance(utctimes, Time):
            utctimes = Time(utctimes)
        return self._astroplan_observer.altaz(utctimes, target = self.target_astroplan)
    
    def risetime(self,
                 utctime : datetime or Time = None ,
                 mode : str = 'nearest',
                 horizon : float = 30,
                 n_grid_points : int = 50) -> Time:
        """
        Calculate the next rise time of the target as seen by the observer.

        Parameters
        ==========
        1. utctime : datetime or Time, optional
            The time to start searching for the next rise time. If not provided, the current time will be used.
        2. mode : str, optional
            The method used to determine the rise time. Possible values are 'next' (the next rise time), 'previous' (the previous rise time), or 'nearest' (the nearest rise time). Default is 'next'.
        3. horizon : float, optional
            The altitude of the horizon, in degrees. Default is 30.

        Returns
        =======
        1. Time
            The rise time of the target as seen by the observer.

        """
        if utctime == None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._astroplan_observer.target_rise_time(utctime, target = self.target_astroplan, which = mode, horizon = horizon*u.deg, n_grid_points = n_grid_points)
    
    def settime(self,
                utctime : datetime or Time or np.array = None,
                mode : str = 'nearest',
                horizon : float = 30,
                n_grid_points : int = 50) -> Time:
        """
        Calculate the time when the target sets below the horizon.

        Parameters
        ==========
        1. utctime : datetime or Time, optional
            The time to use as the reference time for the calculation, by default the current time.
        2. mode : str, optional
            Set to 'nearest', 'next' or 'previous', by default 'nearest'.
        3. horizon : float, optional
            The altitude of the horizon in degrees. Default is 30.

        Returns
        =======
        1. settime : Time
            The time when the target sets below the horizon.
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._astroplan_observer.target_set_time(utctime, self.target_astroplan, which = mode, horizon = horizon*u.deg , n_grid_points = n_grid_points)
    
    def meridiantime(self,
                     utctime : datetime or Time or np.array = None,
                     mode : str = 'nearest',
                     n_grid_points : int = 50) -> Time:
        """
        Calculate the time at which the target passes through the observer's meridian.

        Parameters
        ==========
        1. utctime : datetime or Time, optional
            The time at which to calculate the meridian transit time. If not provided, the current time will be used.
        2. mode : str, optional
            Set to 'nearest', 'next' or 'previous', by default 'nearest'.
            
        Return
        ======
        1. meridiantime : Time
            The time at which the target passes through the observer's meridian.
        """
        
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._astroplan_observer.target_meridian_transit_time(utctime, self.target_astroplan, which = mode, n_grid_points = n_grid_points)
    
    def hourangle(self,
                  utctimes : datetime or Time or np.array = None):
        """
        Calculate the hour angle of the target for a given time(s) in UTC.

        Parameters
        ==========
        1. utctimes : datetime or Time, optional
            The time(s) to calculate the hour angle of the target for, in UTC. If not provided, the current time will be used. 

        Returns
        =======
        1. hourangle : astropy.coordinates.Angle
            The hour angle of the target(s) at the specified time(s).
        """
        
        if utctimes is None:
            utctimes = Time.now()
        if not isinstance(utctimes, Time):
            utctimes = Time(utctimes)
        if not isinstance(self.target_astroplan, FixedTarget):
            raise ValueError('No target is specified for hourangle')
        return self._astroplan_observer.target_hour_angle(utctimes, self.target_astroplan)
    
    def staralt(self,
                utctime : datetime or Time or np.array = None):
        """
        Creates a plot of the altitude and azimuth of a celestial object.
        
        Parameters
        ==========
        1. utctime : datetime or Time or np.array, optional
            The time(s) for which to calculate the altitude and azimuth of the celestial object. 
            If not provided, the current time is used.
        Returns
        =======
        None
        """
        
        now = Time.now()
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)

        sunrisetime = self._observer.tonight(time = utctime, horizon = 0)[1]
        sunsettime = self._observer.sun_settime(sunrisetime, mode = 'previous', horizon= 0)
        astro_sunrisetime = self._observer.sun_risetime(sunrisetime, mode = 'previous', horizon= -18)
        astro_sunsettime = self._observer.sun_settime(sunrisetime, mode = 'previous', horizon= -18)        
        
        #astro_sunrisetime = self._observer.tonight(time = utctime)[1]
        #astro_sunsettime = self._observer.sun_settime(astro_sunrisetime, mode = 'nearest')

        #sunsettime = self._observer.sun_settime(utctime, horizon = 0)
        #sunrisetime = self._observer.sun_risetime(sunsettime, horizon = 0, mode = 'next')
        time_range_start, time_range_end = sunsettime.datetime - datetime.timedelta(hours = 2), sunrisetime.datetime + datetime.timedelta(hours = 2)
        time_axis = np.arange(time_range_start, time_range_end, datetime.timedelta(minutes = 5))
        moon_altaz = self._observer.moon_altaz(time_axis)
        sun_altaz = self._observer.sun_altaz(time_axis)
        target_altaz = self.altaz(time_axis)
        plt.figure(dpi = 300, figsize = (10, 4))
        plt.title(self.name)
        if (now.datetime < time_range_end + datetime.timedelta(hours = 3)) & (now.datetime > time_range_start - datetime.timedelta(hours = 3)):
            plt.axvline(now.datetime, linestyle = '--', c='r', label = 'Now')
        plt.axvline(utctime.datetime, linestyle = '--', c='b', label = 'Observation')
        plt.scatter(moon_altaz.obstime.datetime, moon_altaz.alt.value, c = moon_altaz.az.value, cmap = 'viridis', s = 10, marker = '.', label ='Moon')
        plt.scatter(sun_altaz.obstime.datetime, sun_altaz.alt.value, c = 'k', cmap = 'viridis', s = 15, marker = '.', label = 'Sun')
        plt.scatter(target_altaz.obstime.datetime, target_altaz.alt.value, c = target_altaz.az.value, cmap = 'viridis', s = 30, marker = '*', label = 'Target')
        plt.fill_betweenx([10,90], astro_sunsettime.datetime, astro_sunrisetime.datetime, alpha = 0.1)
        plt.fill_betweenx([10,90], sunsettime.datetime, sunrisetime.datetime, alpha = 0.1)
        plt.axvline(x=astro_sunrisetime.datetime, linestyle = '-', c='k', linewidth = 0.5)
        plt.axvline(x=astro_sunsettime.datetime, linestyle = '-', c='k', linewidth = 0.5)
        plt.axvline(x=sunrisetime.datetime, linestyle = '--', c='k', linewidth = 0.5)
        plt.axvline(x=sunsettime.datetime, linestyle = '--', c='k', linewidth = 0.5)
        #plt.text(astro_sunsettime.datetime, 92, 'Twilight', fontsize = 10)
        plt.text(sunsettime.datetime-datetime.timedelta(minutes=00), 92, 'S.set', fontsize = 10)
        plt.text(sunrisetime.datetime-datetime.timedelta(minutes=00), 92, 'S.rise', fontsize = 10)
        plt.xlim(time_range_start - datetime.timedelta(hours = 1), time_range_end + datetime.timedelta(hours = 1))
        plt.ylim(10, 90)
        plt.legend(loc = 1)
        plt.xlabel('UTC [mm-dd hh]')
        plt.ylabel('Altitude [deg]')
        plt.grid()
        plt.colorbar(label = 'Azimuth [deg]')
        plt.xticks(rotation = 45)
        
    def _get_coordinate_radec(self,
                              ra,
                              dec,
                              frame : str = 'icrs') -> SkyCoord:
        return SkyCoord(ra = ra, dec = dec, frame = frame, unit = (u.deg, u.deg))

    def _get_coordinate_altaz(self,
                              alt,
                              az) -> SkyCoord:
        return SkyCoord(alt = alt, az = az, frame = 'altaz', unit = u.deg)
        
    def _get_target(self,
                    coord,
                    target_name : str = '') -> FixedTarget:
        return FixedTarget(coord = coord, name = target_name)
    
    def _get_constraints(self,
                        TARGET_MINALT : float = None,
                        TARGET_MAXALT : float = None,
                        TARGET_MAX_SUNALT : float = None,
                        TARGET_MOONSEP : float = None,
                        TARGET_MAXAIRMASS : float = None,
                        **kwargs) -> list:
        constraint_all = []
        if (TARGET_MINALT != None) & (TARGET_MAXALT != None):
            constraint_altitude = AltitudeConstraint(min = TARGET_MINALT * u.deg, max = TARGET_MAXALT * u.deg, boolean_constraint = True)
            constraint_all.append(constraint_altitude)
        if (TARGET_MAX_SUNALT != None):
            constraint_atnight = AtNightConstraint(max_solar_altitude = TARGET_MAX_SUNALT * u.deg)
            constraint_all.append(constraint_atnight)
        if TARGET_MOONSEP != None:
            constraint_gallatitude = MoonSeparationConstraint(min = TARGET_MOONSEP * u.deg, max = None)
            constraint_all.append(constraint_gallatitude)
        if TARGET_MAXAIRMASS != None:
            constraint_gallatitude = AirmassConstraint(min = 1, max = TARGET_MAXAIRMASS)
            constraint_all.append(constraint_gallatitude)
        return constraint_all
    
    

# %%

if __name__ == '__main__':
    for name_telescope in ['7DT_01']:
        observer = mainObserver(name_telescope = name_telescope, name_project = 'IMSNG')
        A = mainTarget(name_telescope = name_telescope, name_project = 'IMSNG', observer = observer, target_ra =90, target_dec = -20, target_name = f'ShimTarget[{name_telescope}]')
        A.staralt()
# %%
