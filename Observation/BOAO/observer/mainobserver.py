#%%
import datetime
import pytz
import numpy as np
import matplotlib.pyplot as plt
import inspect
import os
from tqdm import tqdm

from astropy.table import vstack
from astropy.table import Table
from astroplan import AltitudeConstraint, MoonSeparationConstraint
from astroplan import FixedTarget, is_observable
from astroplan import Observer
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, get_sun, get_moon
import astropy.units as u
from astropy.time import Time
from astropy.io import ascii
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import plotly.express as px
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
                 TARGET_MINALT : float = 30,
                 TARGET_MAXALT : float = 90,
                 TARGET_MOONSEP : float = 40,
                 TARGET_PATH : str = '/Users/hhchoi1022/code/hhpy/Observation/BOAO/observer/ObservableTargets.txt',
                 REFSTAR_PATH : str = '/Users/hhchoi1022/code/hhpy/Observation/BOAO/observer/ALLBRICQS_241126_reference.csv',
                 calculate_all_targets : bool = False
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
        self._set_constraints(TARGET_MINALT = TARGET_MINALT, TARGET_MAXALT = TARGET_MAXALT, TARGET_MOONSEP = TARGET_MOONSEP)
        
        self.targets = None
        self.refstars = None
        self._register_target(TARGET_PATH)
        self._register_reference_star(REFSTAR_PATH)
        if calculate_all_targets:
            all_targets = Table()
            if self.targets:
                all_targets = vstack([all_targets, self.targets])
            if self.refstars:
                all_targets = vstack([all_targets, self.refstars])
            
            for target in tqdm(all_targets):
                self._load_altaz_target(target = target, utctime = Time.now())

    def __repr__(self):
        methods = [f'mainObserver.{name}()\n' for name, method in inspect.getmembers(
            mainObserver, predicate=inspect.isfunction) if not name.startswith('_')]
        txt = '[Methods]\n'+''.join(methods)
        return txt
    
    ############ Time ############
    def localtime(self, 
                  utctime : datetime = None):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        localtime : datetime = local utctime for the given site
        """
        if utctime is None:
            utctime = datetime.datetime.utcnow()
        if isinstance(utctime, Time):
            utctime = utctime.datetime
        localtime = pytz.utc.localize(utctime).astimezone(self._timezone)
        return localtime
    
    def siderialtime(self,
                     utctime : datetime.datetime or Time = None,
                     mode : str = 'mean'): # Or apparent 
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        mode : str = mean or apparent, mode to calculate sidereal utctime(default = mean)
        
        Return
        ======
        siderialtime : Longitude = siderialtime utctime for the given site
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._observer.local_sidereal_time(utctime, kind = mode)
        
    def is_night(self,
                 utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        is_night : bool = True if the utctime is night
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._observer.is_night(utctime, horizon = -18*u.deg)
    
    def tonight(self,
                utctime : datetime.datetime or Time = None,
                horizon = -18):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time to calculate tonight(default = Now)
        horizon : float = Degrees above/below actual horizon to use for calculating rise/set times
        
        Return
        ======
        tonight : tuple = tuple[0] for night start utctime, tuple[1] for night end utctime
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        sunrise_civil = self._observer.tonight(utctime, horizon = 0  * u.deg)[1]
        sunrise_night = self.sun_risetime(sunrise_civil, mode = 'previous', horizon= horizon)
        sunset_night = self.sun_settime(sunrise_civil, mode = 'previous', horizon= horizon)        
        return sunset_night, sunrise_night

    ############ Target ############
    
    def to_altaz(self,
                 ra : float or str or list,
                 dec : float or str or list,
                 utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to convert AltAz
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        coordinate : SkyCoord = AltAz coordinate of the target
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        coord = self._observer.altaz(utctime, target = self._get_skycoord(ra, dec))
        return coord.altaz
    
    def to_radec(self,
                 alt : float,
                 az : float,
                 utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(AltAz) coordinates to convert RADec
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        coordinate : SkyCoord = RADec coordinate of the target
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        coord = SkyCoord(alt = alt, az = az, frame ='altaz', location = self._earthlocation, obstime = utctime, unit = 'deg')
        return coord.icrs
    
    def risetime(self,
                 ra : float or str or list,
                 dec : float or str or list,
                 utctime : datetime.datetime or Time = None,
                 mode : str = 'next',
                 horizon : float = 30):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate hourangle
        utctime : datetime.datetime or Time = Time(default = Now)
        mode : {‘next’, ‘previous’, ‘nearest’}(default = next)
        horizon : float = horizon angle in degree(default = -18)
        
        Return
        ======
        risetime : Time = Target risetime
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        if not isinstance(self._get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified')
        rise_time = self._observer.target_rise_time(utctime, self._get_skycoord(ra, dec), which = mode, horizon = horizon*u.deg)
        return rise_time
    
    def settime(self,
                ra : float or str or list,
                dec : float or str or list,
                utctime : datetime.datetime or Time = None,
                mode : str = 'nearest',
                horizon : float = 30):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate hourangle
        utctime : datetime.datetime or Time = Time(default = Now)
        mode : {‘next’, ‘previous’, ‘nearest’}(default = nearest)
        horizon : float = horizon angle in degree(default = -18)
        
        Return
        ======
        settime : Time = Target settime
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        if not isinstance(self._get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified')
        set_time = self._observer.target_set_time(utctime, self._get_skycoord(ra, dec), which = mode, horizon = horizon*u.deg)
        return set_time
    
    def meridiantime(self,
                     ra : float or str or list,
                     dec : float or str or list,         
                     utctime : datetime.datetime or Time = None,
                     mode : str = 'nearest'):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate hourangle
        utctime : datetime.datetime or Time = Time(default = Now)
        mode : {‘next’, ‘previous’, ‘nearest’}(default = nearest)
        
        Return
        ======
        meridiantime : Time = Target meridian transit utctime
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        if not isinstance(self._get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified')
        meridian_time = self._observer.target_meridian_transit_time(utctime, self._get_skycoord(ra, dec), which = mode)
        return meridian_time
    
    def hourangle(self,
                  ra : float or str or list,
                  dec : float or str or list,   
                  utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate hourangle
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        hourangle : Longitude = Hourangle of the target, + For west(pass the meridian)
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        if not isinstance(self._get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified for hourangle')
        return self._observer.target_hour_angle(utctime, self._get_skycoord(ra, dec))                
    
    def parallactic_angle(self,
                          ra : float or str or list,
                          dec : float or str or list,   
                          utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        radec : SkyCoord = Target(RADec) coordinates to calculate parallactic angle
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        parallangle : Angle = parallactic_angle of the target
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        if not isinstance(self._get_skycoord(ra, dec), SkyCoord):
            raise ValueError('No target is specified for hourangle')
        paralangle = self._observer.parallactic_angle(time = utctime, target = self._get_skycoord(ra, dec)).deg
        #if paralangle < 0:
        #    paralangle += 360
        return paralangle

    def is_observable_tonight(self,
                              ra : float, 
                              dec : float,
                              utctime : datetime.datetime or Time or np.array = None,
                              time_grid_resolution : float = 0.1* u.hour) -> bool:
        """
        Determines whether the target is observable at the specified utctime or at the current utctime.

        Parameters
        ----------
        utctime : datetime.datetime or Time, optional
            The utctime at which to check observability. Defaults to the current utctime.
            
        Returns
        -------
        bool
            True if the target is observable, False otherwise.
        """
        tonight = self.tonight(utctime)
        starttime = tonight[0]
        endtime = tonight[1]
        time_range = [starttime, endtime]
        target = FixedTarget(coord = self._get_skycoord(ra, dec))
        return is_observable(constraints = self.constraints, observer = self._observer, targets = target, time_range = time_range, time_grid_resolution = time_grid_resolution)

    ############ Sun ############
    
    def sun_radec(self,
                  utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        sunradec : SkyCoord = RADec coordinate of the Sun
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return get_sun(utctime)
    
    def sun_altaz(self,
                  utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        sunaltaz : SkyCoord = AltAz coordinate of the Sun
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._observer.sun_altaz(utctime)
    
    def sun_risetime(self,
                     utctime : datetime.datetime or Time = None,
                     mode = 'nearest',
                     horizon = -18):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        mode : {‘next’, ‘previous’, ‘nearest’}
        
        
        Return
        ======
        sunrisetime : Time = Sunrise utctime 
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._observer.sun_rise_time(utctime, which = mode, horizon = horizon * u.deg)
    
    def sun_settime(self,
                    utctime : datetime.datetime or Time = None,
                    mode = 'nearest',
                    horizon = -18):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        mode : str = mode {‘next’, ‘previous’, ‘nearest’}
        
        
        Return
        ======
        sunsettime : Time = Sunset utctime 
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._observer.sun_set_time(utctime, which = mode, horizon = horizon * u.deg)
    
    ############ Moon ############
    
    def moon_radec(self,
                   utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        moonradec : SkyCoord = RADec coordinate of the Moon
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return get_moon(utctime)
    
    def moon_altaz(self,
                   utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        moonaltaz : SkyCoord = AltAz coordinate of the Moon
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
            
        return self._observer.moon_altaz(utctime)

    def moon_risetime(self,
                      utctime : datetime.datetime or Time = None,
                      mode = 'nearest',
                      horizon = -18):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        mode : str = mode {‘next’, ‘previous’, ‘nearest’}(default = nearest)
        horizon : float = horizon angle in degree(default = -18)
        
        Return
        ======
        risetime : Time = Moon rise utctime above the horizon
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._observer.moon_rise_time(utctime, which = mode, horizon = horizon * u.deg)
    
    def moon_settime(self,
                     utctime : datetime.datetime or Time = None,
                     mode = 'nearest',
                     horizon = -18):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        mode : str = mode {‘next’, ‘previous’, ‘nearest’}(default = nearest)
        horizon : float = horizon angle in degree(default = -18)
        
        
        Return
        ======
        settime : Time = Moon Set utctime above the horizon
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._observer.moon_set_time(utctime, which = mode, horizon = horizon * u.deg)
    
    def moon_phase(self,
                   utctime : datetime.datetime or Time = None):
        """
        Parameters
        ==========
        utctime : datetime.datetime or Time = Time(default = Now)
        
        Return
        ======
        phase : float = 0 is new, 1 is full
        """
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        return self._observer.moon_illumination(utctime)
    
    ############ Visualization ############
    
    def staralt(self,
                ra : float,
                dec : float,
                objname : str = None,
                utctime : datetime.datetime or Time or np.array = None):
        """
        Creates a plot of the altitude and azimuth of a celestial object.
        
        Parameters
        ----------
        utctime : datetime.datetime or Time or np.array, optional
            The utctime(s) for which to calculate the altitude and azimuth of the celestial object. 
            If not provided, the current utctime is used.
        """
        now = Time.now()
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        astro_sunsettime  = self.sun_settime(utctime, horizon = -18)
        astro_sunrisetime = self.sun_risetime(astro_sunsettime, horizon = -18, mode = 'next')
        sunsettime = self.sun_settime(utctime, horizon = 0)
        sunrisetime = self.sun_risetime(sunsettime, horizon = 0, mode = 'next')
        time_range_start, time_range_end = sunsettime.datetime - datetime.timedelta(hours = 2), sunrisetime.datetime + datetime.timedelta(hours = 2)
        time_axis = np.arange(time_range_start, time_range_end, datetime.timedelta(minutes = 5))
        time_axis_local = np.arange(self.localtime(time_range_start), self.localtime(time_range_end), datetime.timedelta(minutes = 5))
        moon_altaz = self.moon_altaz(time_axis)
        sun_altaz = self.sun_altaz(time_axis)
        target_altaz = self.to_altaz(ra = ra, dec =dec, utctime = time_axis)
        plt.figure(dpi=300, figsize=(10, 4))
        plt.title('Visibility of %s' % objname)

        # Plotting 'Now' line
        if (now.datetime < time_range_end + datetime.timedelta(hours=3)) & (now.datetime > time_range_start - datetime.timedelta(hours=3)):
            plt.axvline(self.localtime(now.datetime), linestyle='--', c='r', label='Now')

        # Plotting Moon, Sun, and Target data
        plt.scatter(time_axis_local, moon_altaz.alt.value, c=moon_altaz.az.value, cmap='viridis', s=10, marker='.', label='Moon')
        plt.scatter(time_axis_local, sun_altaz.alt.value, c='k', cmap='viridis', s=15, marker='.', label='Sun')
        plt.scatter(time_axis_local, target_altaz.alt.value, c=target_altaz.az.value, cmap='viridis', s=30, marker='*', label='Target')

        # Twilight and sunset/sunrise regions
        plt.fill_betweenx([10, 90], self.localtime(astro_sunsettime.datetime), self.localtime(astro_sunrisetime.datetime), alpha=0.1)
        plt.fill_betweenx([10, 90], self.localtime(sunsettime.datetime), self.localtime(sunrisetime.datetime), alpha=0.1)

        # Adding vertical lines for twilight and sunset/sunrise
        plt.axvline(x=self.localtime(astro_sunrisetime.datetime), linestyle='-', c='k', linewidth=0.5)
        plt.axvline(x=self.localtime(astro_sunsettime.datetime), linestyle='-', c='k', linewidth=0.5)
        plt.axvline(x=self.localtime(sunrisetime.datetime), linestyle='--', c='k', linewidth=0.5)
        plt.axvline(x=self.localtime(sunsettime.datetime), linestyle='--', c='k', linewidth=0.5)

        # Adding labels for twilight, sunset, and sunrise
        plt.text(self.localtime(astro_sunsettime.datetime) - datetime.timedelta(minutes=0), 92, 'Twilight', fontsize=10)
        plt.text(self.localtime(sunsettime.datetime) - datetime.timedelta(minutes=0), 92, 'S.set', fontsize=10)
        plt.text(self.localtime(sunrisetime.datetime) - datetime.timedelta(minutes=0), 92, 'S.rise', fontsize=10)

        # Setting x and y limits
        plt.xlim(self.localtime(time_range_start - datetime.timedelta(hours=1)), self.localtime(time_range_end + datetime.timedelta(hours=1)))
        plt.ylim(10, 90)

        plt.legend(loc=1)
        plt.xlabel('Local Time [mm-dd hh]')
        plt.ylabel('Altitude [deg]')
        plt.grid()
        plt.colorbar(label='Azimuth [deg]')
        plt.xticks(rotation=45)
        plt.show()
            
    def visualize_altaz(self, utctime: datetime.datetime or Time = None, update_time: float = 10, visualize_track: bool = True):
        """
        Visualize altitude and azimuth in a polar plot using Dash with dynamic updates.
        - Targets and reference stars are displayed with different colors.
        - Hovering over points shows detailed information.
        - The plot updates dynamically every `update_time` seconds.
        """
        # Create a Dash app
        app = dash.Dash(__name__)

        # Helper function to convert Alt-Az to polar coordinates
        def _convert_altaz_to_polar(alt, az):
            r = 90 - alt  # Altitude -> Distance from center (0 = zenith)
            theta = az  # Azimuth in degrees
            return theta, r        
            
        def _generate_figure():
            fig = go.Figure()

            # Define the range for the azimuth (0 to 360 degrees)
            theta = np.linspace(0, 360, 360)

            # Define the radius for the outer boundary of the plot (e.g., r = 90)
            r_outer = np.full_like(theta, 90)

            # Define the radius for the inner boundary of the filled region (r = 60)
            r_inner = np.full_like(theta, 60)

            # Add a filled region for r > 60
            fig.add_trace(
                go.Scatterpolar(
                    theta=np.concatenate([theta, theta[::-1]]),  # Combine for the fill
                    r=np.concatenate([r_outer, r_inner[::-1]]),  # Outer and inner radius
                    fill='toself',  # Fill between the two boundaries
                    fillcolor='rgba(255, 0, 0, 0.1)',  # Semi-transparent red
                    line=dict(color='red'),  # Outline color
                    name='Region Alt < 30'
                )
            )

            # Update targets
            if self.targets is not None:
                targets = self.targets.copy()
                altaz_targets = self.to_altaz(targets['ra'], targets['dec'], utctime=utctime)
                targets['alt'] = np.round(altaz_targets.alt.deg, 2)
                targets['az'] = np.round(altaz_targets.az.deg, 2)
                targets['parallactic_angle'] = np.round(self.parallactic_angle(targets['ra'], targets['dec'], utctime), 2)

                theta_targets, r_targets = _convert_altaz_to_polar(targets['alt'], targets['az'])

                # Add target markers
                txt_for_each_target = [
                    ''.join([f'{key}: {value}<br>' for key, value in dict(target).items()])
                    for target in targets
                ]

                fig.add_trace(
                    go.Scatterpolar(
                        theta=theta_targets,
                        r=r_targets,
                        mode="markers",
                        marker=dict(
                            size=10,
                            color=["blue" if status == "OBSERVED" else "black" for status in targets["status"]],
                            symbol=['star' if priority < 1.1 else 'circle' for priority in targets["priority"]]
                        ),
                        name="Targets",
                        text=txt_for_each_target,
                        hoverinfo="text",
                    )
                )

                # Add target tracks using plotly.express.line_polar
                if visualize_track:
                    targets_above_horizon = targets[targets['alt'] > 0]
                    for target in targets_above_horizon:
                        target_altaz_tbl = self._load_altaz_target(target=target, utctime=utctime)
                        alt = target_altaz_tbl['alt']
                        az = target_altaz_tbl['az']
                        theta_track, r_track = _convert_altaz_to_polar(alt, az)
                        track_tbl = Table({'theta': theta_track, 'r': r_track}) 
                        
                        # Use plotly.express.line_polar for drawing the tracks
                        import pandas as pd
                        track_df = track_tbl.to_pandas()
                        track_fig = px.line_polar(
                            track_df,
                            r='r',
                            theta='theta',
                            line_close=False,
                        )
                        track_fig.update_traces(line=dict(width=0.1, color='grey'))
                        fig.add_trace(track_fig.data[0])

            # Update reference stars
            if self.refstars is not None:
                refstars = self.refstars.copy()
                altaz_refstars = self.to_altaz(refstars["ra"], refstars["dec"], utctime=utctime)
                refstars["alt"] = np.round(altaz_refstars.alt.deg, 2)
                refstars["az"] = np.round(altaz_refstars.az.deg, 2)

                theta_refstars, r_refstars = _convert_altaz_to_polar(refstars["alt"], refstars["az"])

                # Add reference star markers
                txt_for_each_refstar = [
                    ''.join([f'{key}: {value}<br>' for key, value in dict(refstar).items()])
                    for refstar in refstars
                ]

                fig.add_trace(
                    go.Scatterpolar(
                        theta=theta_refstars,
                        r=r_refstars,
                        mode="markers",
                        marker=dict(size=10, color="red"),
                        name="Reference Stars",
                        text=txt_for_each_refstar,
                        hoverinfo="text",
                    )
                )
                
                # Add restarts tracks using plotly.express.line_polar
                if visualize_track:
                    refstars_above_horizon = refstars[refstars['alt'] > 0]
                    for refstar in refstars_above_horizon:
                        refstar_altaz_tbl = self._load_altaz_target(target=refstar, utctime=utctime)
                        alt = refstar_altaz_tbl['alt']
                        az = refstar_altaz_tbl['az']
                        theta_track, r_track = _convert_altaz_to_polar(alt, az)
                        track_tbl = Table({'theta': theta_track, 'r': r_track}) 
                        
                        # Use plotly.express.line_polar for drawing the tracks
                        import pandas as pd
                        track_df = track_tbl.to_pandas()
                        track_fig = px.line_polar(
                            track_df,
                            r='r',
                            theta='theta',
                            line_close=False,
                        )
                        track_fig.update_traces(line=dict(width=0.1, color='grey'))
                        fig.add_trace(track_fig.data[0])


            # Update layout
            fig.update_layout(
                polar=dict(
                    angularaxis=dict(
                        direction="clockwise",
                        tickmode="array",
                        tickvals=[0, 270, 180, 90],  # Azimuth angles for "N", "E", "S", "W"
                        ticktext=["N", "W", "S", "E"],  # Cardinal directions
                    ),
                    radialaxis=dict(
                        range=[0, 90],  # Altitude range: 0 (horizon) to 90 (zenith)
                        tickmode="array",
                        tickvals=[90 - val for val in [15, 30, 45, 60, 75, 90]],  # Custom altitude ticks
                        ticktext=[f"{val}°" for val in [15, 30, 45, 60, 75, 90]],  # Add degree symbols
                    ),
                ),
                showlegend=True,
            )

            return fig


        # Layout of the Dash app
        app.layout = html.Div([
            html.Div(id='current-time', style={'textAlign': 'center', 'fontSize': 18, 'marginBottom': '10px'}),
            dcc.Graph(id='polar-plot', style={'height': '90vh'}),
            dcc.Interval(id='interval-component', interval=update_time * 1000, n_intervals=0)  # Interval in milliseconds
        ])

        # Callback to update the figure
        @app.callback(
            [Output('polar-plot', 'figure'), Output('current-time', 'children')],
            [Input('interval-component', 'n_intervals')]
        )
        def update_plot(n):
            # Update time
            current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            return _generate_figure(), f"Updated Time: {current_time}"

        # Run the Dash app
        app.run_server(debug=True, use_reloader=False, port = 8049)
        print('Dash app is running at http://127.0.0.1:8049/')

    ############ Additional functions ############
    # Find optimal reference star for the given coordinates
    def get_reference_star(self,
                           ra : float or str,
                           dec : float or str,
                           utctime : datetime.datetime or Time = None,
                           idx_from_smallest : int = 0,
                           sort_by : str = 'separation',
                           altitude_cut : float = 30,
                           consider_meridian_flip: bool = True,
                           meridian_flip_angle_tolerance: float = 15):
        if self.refstars is None:
            print('[WARNING] No reference star table is registered')
            return None
            
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
            
        target_coord = self._get_skycoord(ra, dec)
        target_altaz = self.to_altaz(target_coord.ra.deg, target_coord.dec.deg, utctime = utctime)
        refstar_tbl = self.refstars.copy()
        refstar_coord = self._get_skycoord(refstar_tbl['ra'], refstar_tbl['dec'])
        refstar_altaz = self.to_altaz(refstar_coord.ra.deg, refstar_coord.dec.deg, utctime = utctime)
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
                filter_refstar  = refstar_tbl['az'] > 180
            else:
                filter_refstar  = refstar_tbl['az'] < 180 - meridian_flip_angle_tolerance
            refstar_tbl = refstar_tbl[filter_refstar]
            if len(refstar_tbl) == 0:
                raise ValueError('No reference star is found due to meridian flip criteria')
            
        refstar = refstar_tbl[idx_from_smallest]
        print('Target Alt = %.1f, Az = %.1f'%(target_altaz.alt.deg, target_altaz.az.deg))
        print('Ref Alt = %.1f, Az = %.1f'%(refstar['alt'], refstar['az']))
        print('Ref Separation = %.1f'%refstar['separation'])
        return refstar, refstar_tbl
        
    # update target table
    def update_target_table(self, objname: str, value : str, column_key: str = 'status'):
        """
        Updates the target table with a given value for a specific objname and column_key.

        Parameters
        ----------
        objname : str
            The name of the target object to update.
        column_key : str
            The column key in the target table to update.
        value : Any
            The new value to assign to the specified column for the specified objname.
        """
        # Ensure the column exists
        if column_key not in self.targets.colnames:
            raise ValueError(f"Column '{column_key}' does not exist in the target table.")

        # Get the index of the row corresponding to the objname
        index = np.where(self.targets['objname'] == objname)[0]

        if len(index) == 0:
            raise ValueError(f"Object '{objname}' not found in the target table.")

        # Update the value explicitly using the index
        self.targets[column_key][index[0]] = value
        print(f"Updated '{objname}' in column '{column_key}' with value '{value}'.")

    # update reference star table
    def update_reference_star_table(self, objname: str, value : str, column_key: str = 'status'):
        """
        Updates the reference star table with a given value for a specific objname and column_key.

        Parameters
        ----------
        objname : str
            The name of the reference star object to update.
        column_key : str
            The column key in the reference star table to update.
        value : Any
            The new value to assign to the specified column for the specified objname.
        """
        # Ensure the column exists
        if column_key not in self.refstars.colnames:
            raise ValueError(f"Column '{column_key}' does not exist in the reference star table.")

        # Get the index of the row corresponding to the objname
        index = np.where(self.refstars['Name'] == objname)[0]

        if len(index) == 0:
            raise ValueError(f"Object '{objname}' not found in the reference star table.")

        # Update the value explicitly using the index
        self.refstars[column_key][index[0]] = value
        print(f"Updated '{objname}' in column '{column_key}' with value '{value}'.")
        
    def _register_reference_star(self,
                                 refstar_path : str = './ALLBRICQS_241126_reference.csv'):
        try:
            refstar_tbl = ascii.read(refstar_path)
            coord = self._get_skycoord(refstar_tbl['ra'], refstar_tbl['dec'])
            refstar_tbl['ra'] = np.round(coord.ra.deg,4)
            refstar_tbl['dec'] = np.round(coord.dec.deg,4)
            coords = self._get_skycoord(refstar_tbl['ra'], refstar_tbl['dec'])
            refstar_tbl['ra_hms'] = coords.ra.to_string(unit=u.hourangle, sep=' ', precision=2, pad=True)
            refstar_tbl['dec_dms'] = coords.dec.to_string(unit=u.deg, sep=' ', precision=2, pad=True, alwayssign = True)
            refstar_tbl['note'] = ''
            self.refstars = refstar_tbl
        except:
            raise print('[WARNING] Reference star table cannot be loaded')

    def _register_target(self,
                         target_path : str = './ObservableTargets.txt',
                         ):
        try:
            target_tbl = ascii.read(target_path)
            coords = self._get_skycoord(target_tbl['ra'], target_tbl['dec'])
            target_tbl['ra_hms'] = coords.ra.to_string(unit=u.hourangle, sep=' ', precision=2, pad=True)
            target_tbl['dec_dms'] = coords.dec.to_string(unit=u.deg, sep=' ', precision=2, pad=True, alwayssign = True)
            target_tbl['exptime'] = 0
            target_tbl['status'] = 'NOT_OBSERVED'
            target_tbl['note'] = ''
            target_tbl['is_observable_tonight'] = self.is_observable_tonight(ra = target_tbl['ra'], dec = target_tbl['dec'], utctime = Time.now())
            self.targets = target_tbl
        except:
            raise print('[WARNING] Target table cannot be loaded')
    
    def _get_skycoord(self,
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
    
    def _set_constraints(self,
                         TARGET_MINALT : float = 30,
                         TARGET_MAXALT : float = 90,
                         TARGET_MOONSEP : float = 40,
                         **kwargs) -> list:
        constraint_all = []
        if (TARGET_MINALT != None) & (TARGET_MAXALT != None):
            constraint_altitude = AltitudeConstraint(min = TARGET_MINALT * u.deg, max = TARGET_MAXALT * u.deg, boolean_constraint = True)
            constraint_all.append(constraint_altitude)
        if TARGET_MOONSEP != None:
            constraint_gallatitude = MoonSeparationConstraint(min = TARGET_MOONSEP * u.deg, max = None)
            constraint_all.append(constraint_gallatitude)
        self.constraints = constraint_all
        
    def _load_altaz_target(self, target, utctime : datetime.datetime or Time = None):
        target = Table(target)
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
        def _calculate_and_save_track_target(target, utctime):
            sunset, sunrise = self.tonight(utctime)
            time_range = np.arange(sunset.datetime, sunrise.datetime, datetime.timedelta(minutes = 10))
            altaz_track = self.to_altaz(target['ra'], target['dec'], utctime = time_range)
            track_altaz_tbl = Table()
            track_altaz_tbl['utctime'] = time_range
            track_altaz_tbl['alt'] = altaz_track.alt.deg
            track_altaz_tbl['az'] = altaz_track.az.deg
            track_altaz_tbl.write(os.path.join('./_data', 'altaz', Time(utctime).datetime.strftime('%Y%m%d '), f'trackinfo_{target["objname"][0]}.csv'))
        # Set trac
        track_path = os.path.join('./_data', 'altaz', Time(utctime).datetime.strftime('%Y%m%d '), f'trackinfo_{target["objname"][0]}.csv')
        
        # if the track folder does not exist, make the folder
        if not os.path.exists(os.path.dirname(track_path)):
            os.makedirs(os.path.dirname(track_path), exist_ok=True)
        
        # if the track file does not exist, calculate and save the track
        if not os.path.exists(track_path):
            _calculate_and_save_track_target(target, utctime)
        
        # load the track file
        track_info = ascii.read(track_path)
        return track_info
        
# %% fractical use
if __name__ == '__main__':
    observer = mainObserver(calculate_all_targets= True)
#%%
if __name__ == '__main__':
    utctime = Time.now()
    localtime = observer.localtime(utctime)
    sunset, sunrise = observer.tonight(utctime)
    sunset_ltc = observer.localtime(sunset)
    sunrise_ltc = observer.localtime(sunrise)
    moon_phase = observer.moon_phase(utctime)
    moon_altaz = observer.moon_altaz(utctime)
    print('Currnet Time: ', localtime)
    print('Sunset Time: ', sunset_ltc.strftime('%Y-%m-%d %H:%M:%S'))
    print('Sunrise Time: ', sunrise_ltc.strftime('%Y-%m-%d %H:%M:%S'))
    print('Moon Phase: ', moon_phase)
    print('Moon Alt/Az: ', np.round(moon_altaz.alt.deg,2), np.round(moon_altaz.az.deg,2))
    
    target_name = 'J0218+4534'
    ra, dec = observer.targets[observer.targets['objname'] == target_name]['ra'], observer.targets[observer.targets['objname'] == target_name]['dec']
    risetime = observer.risetime(ra, dec, utctime = utctime, mode = 'nearest', horizon = 35)
    meridian_time = observer.meridiantime(ra, dec, utctime = risetime, mode ='next')
    settime = observer.settime(ra, dec, utctime = risetime, mode ='next', horizon = 35)
    risetime_ltc = observer.localtime(risetime)
    meridian_time_ltc = observer.localtime(meridian_time)
    settime_ltc = observer.localtime(settime)
    parallel_angle = observer.parallactic_angle(ra, dec, utctime = utctime)
    hourangle = observer.hourangle(ra, dec, utctime = utctime)
    is_observable_tonight = observer.is_observable_tonight(ra, dec, utctime = Time.now())
    refstar, refstar_all = observer.get_reference_star(ra,
                                                       dec,
                                                       utctime= None,
                                                       idx_from_smallest = 0,
                                                       sort_by = 'separation',
                                                       altitude_cut = 30,
                                                       consider_meridian_flip = True,
                                                       meridian_flip_angle_tolerance = 15)
    
    print('Target: ', target_name)
    print('RA%.5f/Dec%.5f'%(ra[0], dec[0]))
    print('Rise Time: ', risetime_ltc.strftime('%Y-%m-%d %H:%M:%S'))
    print('Meridian Time: ', meridian_time_ltc.strftime('%Y-%m-%d %H:%M:%S'))
    print('Set Time: ', settime_ltc.strftime('%Y-%m-%d %H:%M:%S'))
    print('Parallactic Angle: %.2f'%(parallel_angle))
    print('Hourangle: %.2f '%(hourangle[0].value))
    print('Is_Observable_tonight: ', is_observable_tonight[0])
    print('Reference Star: ', refstar)
# %%
observer.staralt(ra, dec, objname = target_name, utctime = utctime)
#%%
observer.update_target_table('J1811+2409', 'OBSERVED', column_key = 'status')
observer.update_target_table('J0425+3454', 'OBSERVED', column_key = 'status')
observer.update_target_table('J0422+2854', 'OBSERVED', column_key = 'status')
observer.update_target_table('J2206+2757', 'OBSERVED', column_key = 'status')
observer.update_target_table('J1918+4222', 'OBSERVED', column_key = 'status')
observer.update_target_table('J2139+2929', 'OBSERVED', column_key = 'status')
observer.update_target_table('J2135+3348', 'OBSERVED', column_key = 'status')
observer.update_target_table('J2225+3952', 'OBSERVED', column_key = 'status')
observer.update_target_table('J0214+3824', 'OBSERVED', column_key = 'status')
observer.update_target_table('J0218+4534', 'OBSERVED', column_key = 'status')
observer.update_target_table('J0319+3309', 'OBSERVED', column_key = 'status')
observer.update_target_table('J1846+8425', 'OBSERVED', column_key = 'status')
observer.update_target_table('J0422+2854', 'OBSERVED', column_key = 'status')
observer.update_target_table('J0425+3454', 'OBSERVED', column_key = 'status')
observer.update_target_table('J1505+8604', 'OBSERVED', column_key = 'status')
observer.update_target_table('J1317+8637', 'OBSERVED', column_key = 'status')
#observer.update_reference_star_table('HZ', 'target_name', column_key = 'note')
# %%
observer.visualize_altaz(visualize_track = True)
#%%

#%
data = ascii.read('./241126_stdstars_eso.csv', delimiter=',')
data.rename_columns(['col1','col2','col3','col4','col5','col6'],['Name','ra','dec','V','type','note'])

coords = observer._get_skycoord(data['ra'], data['dec'])
data['ra_hms'] = coords.ra.to_string(unit=u.hourangle, sep=' ', precision=2, pad=True)
data['dec_dms'] = coords.dec.to_string(unit=u.deg, sep=' ', precision=2, pad=True, alwayssign = True)
# %%
from astropy.table import join
# %%
data
# %%
data_coord = ascii.read('./refstar_1.csv', delimiter = ',')
data_coord.rename_columns(['col1','col2','col3'],['Name','ra','dec'])
data_magnitude = ascii.read('refstar_2.csv', delimiter = ',')
data_magnitude.rename_columns(['col1','col2','col3','col4'],['Name','type','V','B-V'])
# %%
merged_tbl = join(left = data_coord, right = data_magnitude, keys = 'Name', join_type = 'inner')
merged_tbl.remove_columns(['col4','col5_1','col6_1','col7_1','col8_1','col5_2','col6_2','col7_2','col8_2'])
# %%
merged_tbl = merged_tbl[merged_tbl['V'].astype(float) < 14]
merged_tbl.write('ALLBRICQS_241126_reference.csv', format = 'csv', overwrite = True)
# %%
