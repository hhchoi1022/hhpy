#%%
from mainobserver import mainObserver
import datetime
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
#%%
class Staralt():
    
    def __init__(self, 
                 observer : mainObserver,
                 utctime : datetime or Time = None,
                 ):
        # Set the observer
        self._observer = observer
        # If no time is provided, use the current time
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        # Set the night
        self.tonight = self._set_night(utctime = utctime)
        
    def _set_night(self, utctime : datetime or Time):
        """
        Set the night for the given time.
        
        Parameters
        ----------
        utctime : datetime or Time
            The time for which to set the night.
        """
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        class obsnight: 
            def __repr__(self):
                attrs = {name: value.iso if isinstance(value, Time) else value
                         for name, value in self.__dict__.items()}
                max_key_len = max(len(key) for key in attrs.keys())
                attrs_str = '\n'.join([f'{key:{max_key_len}}: {value}' for key, value in attrs.items()])
                return f'{self.__class__.__name__} Attributes:\n{attrs_str}'
        obsnight = obsnight()
        # Celestial information
        obsnight.sunrise_civil = self._observer.sun_risetime(utctime, horizon = 0, mode = 'next')
        obsnight.sunset_civil = self._observer.sun_settime(obsnight.sunrise_civil, mode = 'previous', horizon= 0)        
        obsnight.sunrise_nautical = self._observer.sun_risetime(obsnight.sunrise_civil, mode = 'previous', horizon= -6)
        obsnight.sunset_nautical = self._observer.sun_settime(obsnight.sunrise_civil, mode = 'previous', horizon= -6)
        obsnight.sunrise_astro = self._observer.sun_risetime(obsnight.sunrise_civil, mode = 'previous', horizon= -12)
        obsnight.sunset_astro = self._observer.sun_settime(obsnight.sunrise_civil, mode = 'previous', horizon= -12)
        obsnight.sunrise_night = self._observer.sun_risetime(obsnight.sunrise_civil, mode = 'previous', horizon= -18)
        obsnight.sunset_night = self._observer.sun_settime(obsnight.sunrise_civil, mode = 'previous', horizon= -18) 
        return obsnight
        
    def staralt(self, 
                utctime : datetime or Time = None,
                # Target information
                ra : float or str = None,
                dec : float or str = None,
                objname : str = None,
                target_minalt : float = 30,
                target_minmoonsep : float = 40,
                ):
        """
        Creates a plot of the altitude and azimuth of a celestial object.
        
        Parameters
        ----------
        utctime : datetime or Time or np.array, optional
            The time(s) for which to calculate the altitude and azimuth of the celestial object. 
            If not provided, the current time is used.
        """
        # If no time is provided, use the current time
        now = Time.now()
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        
        tonight = self.tonight
        # Calculate the coordinates of the sky objects
        coord = self._get_skycoord(ra, dec)
        time_range_start, time_range_end = tonight.sunset_astro.datetime - datetime.timedelta(hours = 2), tonight.sunrise_astro.datetime + datetime.timedelta(hours = 2)
        time_axis = np.arange(time_range_start, time_range_end, datetime.timedelta(minutes = 5))
        moon_altaz = self._observer.moon_altaz(time_axis)
        sun_altaz = self._observer.sun_altaz(time_axis)
        target_altaz = self._observer.to_altaz(ra = ra, dec = dec, utctime = time_axis)
        target_moonsep = moon_altaz.separation(target_altaz)
        
        # Plot
        plt.figure(dpi = 300, figsize = (10, 4))
        if objname is None:
            titlename = f'Altitude of the Target'
        else:
            titlename = f'Altitude of the {objname}'
        plt.title(titlename)
        if (now.datetime < time_range_end + datetime.timedelta(hours = 3)) & (now.datetime > time_range_start - datetime.timedelta(hours = 3)):
            plt.axvline(now.datetime, linestyle = '--', c='r', label = 'Now')
        color_target_alt = ['k' if alt > target_minalt else 'r' for alt in target_altaz.alt.value]
        color_target_moonsep = ['k' if sep > target_minmoonsep else 'r' for sep in target_moonsep.value]
        color_target = ['r' if 'r' in (x, y) else 'k' for x, y in zip(color_target_alt, color_target_moonsep)]
        plt.scatter(moon_altaz.obstime.datetime, moon_altaz.alt.value, c = 'b', s = 10, marker = '.', label ='Moon')
        plt.scatter(sun_altaz.obstime.datetime, sun_altaz.alt.value, c = 'r', s = 15, marker = '.', label = 'Sun')
        plt.scatter(target_altaz.obstime.datetime, target_altaz.alt.value, c = color_target, s = 30, marker = '*', label = 'Target')
        plt.fill_betweenx([10,90], tonight.sunset_night.datetime, tonight.sunrise_night.datetime, color = 'k', alpha = 0.3)
        plt.fill_betweenx([10,90], tonight.sunset_civil.datetime, tonight.sunrise_civil.datetime, color = 'k', alpha = 0.1)
        plt.axvline(x=tonight.sunrise_night.datetime, linestyle = '-', c='k', linewidth = 0.5)
        plt.axvline(x=tonight.sunset_night.datetime, linestyle = '-', c='k', linewidth = 0.5)
        plt.fill_between([tonight.sunset_night.datetime, tonight.sunrise_night.datetime], 0, target_minalt, color = 'r', alpha = 0.3)
        plt.text(tonight.sunset_night.datetime, 93, 'Night start', fontsize = 10, ha='center', va='center')
        plt.text(tonight.sunrise_night.datetime, 93, 'Night end', fontsize = 10, ha='center', va='center')
        plt.text((tonight.sunset_night + 0.5 * (tonight.sunrise_night.jd - tonight.sunset_night.jd)).datetime, 20, 'Observation limit', fontsize = 10, ha='center', va='center', fontweight='bold', c = 'darkred')
        plt.text(time_range_start - datetime.timedelta(hours = 0.7), 85, f'Current observation criteria:\n- Altitude > {target_minalt} deg\n- Moon separation > {target_minmoonsep} deg', fontsize = 10, ha='left', va='top', c = 'k', fontweight='bold')
        plt.xlim(time_range_start - datetime.timedelta(hours = 1), time_range_end + datetime.timedelta(hours = 1))
        plt.ylim(10, 90)
        plt.legend(loc = 1)
        plt.xlabel('UTC [mm-dd hh]')
        plt.ylabel('Altitude [deg]')
        plt.grid()
        plt.xticks(rotation = 45)

    def _get_skycoord(self,
                     ra : str or float,
                     dec: str or float):
        """
        Parameters
        ==========
        ra : str | float = Right Ascension, if str, it should be in hms format (e.g., "10:20:30"), if float, it should be in decimal degrees
        dec : str | float = Declination, if str, it should be in dms format (e.g., "+20:30:40"), if float, it should be in decimal degrees
        
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
            raise ValueError("Unsupported RA and Dec format")
        return coord
# %% Example 
if __name__ == '__main__':
    import time
    A = mainObserver()
    S = Staralt(A)
    start = time.time()
    S.staralt(ra = 20.5243, dec = -20.245, objname = 'ABC', target_minalt = 30)
    print(f'Elapsed time: {time.time() - start:.2f} s')
# %%
