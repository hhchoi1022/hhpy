

#%%

from mainconfig import mainConfig
from mainobserver import mainObserver
from maintarget import mainTarget
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table

from astropy.io import ascii
import astropy.units as u

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
import uuid
import os

from astroplan import AltitudeConstraint, AirmassConstraint, MoonSeparationConstraint, AtNightConstraint
from astroplan import observability_table

#%%
class ObsScheduler(mainConfig):
    """
    An observation scheduler class that generates observing schedules based on a target database, time constraints, and the telescope name.
    - History
    (23.04.22) Written by Hyeonho Choi 
    (23.04.23) Change in write_ACPscript_RASA36, write_ACPscript_KCT : change the order to apply #NOSOLVING for all targets
    (23.04.25) ""Version 2.0"", Scheduled table data format changed from Astropy.Table to Pandas.Dataframe
    (23.04.26) Added write_txt function
    (23.04.26) Change the output format of write_txt function
    (23.04.26) Change the structure of the configuration directory depending on the project (ToO, IMSNG)
    (23.04.27) No error anymore when no observable target exists
    (23.04.27 11:10) #Minor change // make_ACPscript_LSGT, fix for the case no target exist in split_table
    (23.05.02) ""Version 3.0"", Running time decreased.
    (23.05.03) Fix for config setting (name_project added)
    (23.06.02) (critical) Fix for the coordinate transformation to string (00:30:30 -> -00:30:30)
    (23.07.10) NINA script converter added by Hongjae Moon
    (23.11.01) Major change in the structure 
        
    ==========
    Parameters
    ==========
    1. target_db : Astropy.table.Table, pandas.Dataframe, Dict
        data containing target information.
        - Must indluded keys
        --------
        obj : object name
        ra : target right ascension, '10:00:00', "10 00 00", "10h 00m 00s", '150.0000'
        dec : target declination, '-30:00:00', "-30 00 00", "-30d 00m 00s", '-30.0000'
        Other keys (filters, exptime, ...) will be set default(from Scheduler_telescope.config) if not provided
        --------
        
    2. date : datetime.datetime, astropy.time.Time, str (default=Time.now())
        The observation date (in UT).
        - Current supported string format
        --------
        str : 20230421, 2023/04/21, 2023.04.21, 230421, 23/04/21, 23.04.21
        astropy.time.Time
        datetime.datetime
        --------
        
    3. project : str (default='ToO')
        The name of the project. For GECKO, just use ToO

    4. name_telescope : str (default='KCT')
        The name of the telescope for which the observation plan is being made.
        - Current supported telescopes
        --------
        (ToO project) name_telescope : 'KCT', 'RASA36', 'CBNUO', 'LSGT', 'SAO', 'LOAO' 'KMTNet_CTIO', 'KMTNet_SAAO', 'KMTNet_SSO'
        (IMSNG project) name_telescope : 'KCT', 'RASA36', 'LSGT', 'CBNUO'
        --------

    5. entire_night : bool (default=False)
        Whether to generate a plan for the entire night or only part of it.
    
    =======
    Methods
    =======
    1. scheduler(**kwargs) -> schedule : object
        Schedule the observing plan with the observable target
        before running scheduler(), self.schedule returns None
        
    2. scorer(obs_tbl : Astropy.table.Table, utctime : datetime or Astropy.time.Time, **kwargs) -> score : np.array
        Return score array for the given targetlist(obs_tbl)
        
    3. show(scheduled_table : astropy.table.Table, **kwargs) -> None
        Visualize the observing plan with the scheduled table 
    
    ==========
    Properties
    ==========
    1. constraint 
        - maxalt, minalt, moon_separation
    2. self.observer : mainObserver class
    3. self.obsnight
        - sunset_astro, sunrise_astro : -18deg sunrise/set time
        - sunset_prepare, sunrise_prepare : -5deg sunrse/set time
        - obs_start, obs_end : observation start/ end time. 
        - midnight : observatio midnight
    4. self.obsinfo
        - is_night : night status of the observing site 
        - moon_phase : moon phase at the site and date
        - observer_info : observing site information
        - sun_radec : RADec of the Sun
        - moon_radec : RADec of the Moon
    5. self.target
        - all : all targetlist input
        - observable : observable targetlist
    6. self.schedule
        Only when scheduler() is already run
        - all: all targetlist with the updated "scheduled" status
        - observable: observable targetlist with the updated "scheduled" status
        - scheduled: scheduled targetlist 
    
    """
    def __init__(self,
                    target_db,
                    date = Time.now(),
                    name_project : str = 'ToO',
                    name_telescope = 'KCT',
                    entire_night : bool = False,
                    ):
        super().__init__(name_telescope = name_telescope, name_project = name_project)
        
        self.name_telescope = name_telescope
        self.name_project = name_project
        self._date = date
        self._entire_night = entire_night
        self._target_db = target_db
        
        self.constraint = self._set_constrints()
        self.observer = mainObserver(name_telescope = name_telescope, name_project = name_project)
        self.obsnight = self._set_obs_night(horizon_prepare=-5, horizon_astro = -15)
        self.obsinfo = self._set_obs_info()
        self.target = self._set_target()
        self.schedule = None

    def scheduler(self,
                  duplicate_when_empty : bool = False,
                  
                  autofocus : bool = True,
                  autofocus_at_init : bool = False,
                  autofocus_interval_in_hour : float = 3,
                  
                  n_target_for_each_timeslot : int = 1,
                  delay_minute : float = 10,
                  **kwargs):
        """
        duplicate_when_empty : bool = False
        
        autofocus : bool = True
        autofocus_at_init : bool = False
        autofocus_interval_in_hour : float = 3
        
        n_target_for_each_timeslot : int = 1
        delay_minute : float = 10
        """
        class schedule: pass
        schedule.obsnight = self.obsnight
        schedule.all = self.target.all.to_pandas()
        if len(self.target.observable) > 0:
            scheduled_tbl = Table()
            observable_tbl = self.target.observable
            if 'obstime_start' in self.target.observable.colnames:
                scheduled_tbl = self.target.observable[self.target.observable['obstime_start'] !='']
                observable_tbl = self.target.observable[self.target.observable['obstime_start'] =='']
                scheduled_tbl['obstime_start'] = Time(scheduled_tbl['obstime_start'])
                if len(scheduled_tbl) > 0:
                    scheduled_tbl = scheduled_tbl.to_pandas()
                else:
                    scheduled_tbl = Table().to_pandas()
                if len(observable_tbl) == 0:
                    observable_tbl = self.target.observable[0:1]
            schedule.observable = observable_tbl.to_pandas()
            schedule.scheduled = schedule.observable.copy()[0:0]
            schedule.unscheduled = schedule.observable[schedule.observable['scheduled'] == False]

            def calculate_tot_exptime(target):
                exptimes = str(target['exptime']).split(',')
                counts = str(target['count']).split(',')
                tot_exptime = 0
                for exptime, count in zip(exptimes, counts):
                    tot_exptime += (float(exptime) * float(count))
                return tot_exptime
            
            def calculate_tot_overheadtime(target,
                                        sepdist_slew : float = 60 #deg
                                        ):
                filters = str(target['filter']).split(',')
                counts = str(target['count']).split(',')
                num_totimg = np.array(counts).astype(float).sum()
                num_filt = len(filters)
                time_readout = self.config['OVERHEAD_TIME_READOUT'] * num_totimg
                time_filtchange = self.config['OVERHEAD_TIME_FILTCHAN'] * num_filt
                time_slew = (sepdist_slew/self.config['OVERHEAD_SPEED_SLEWING']) * 3
                time_overhead = time_readout + time_filtchange + time_slew
                return time_overhead
            
            def insert_schedule(obs_start_time,
                                scheduled_table,
                                target,
                                score = None
                                ):
                output_table = scheduled_table.copy()
                target_tmp = target.copy()
                obs_start_time = Time(obs_start_time)
                tot_exptime = calculate_tot_exptime(target = target)
                end_time = obs_start_time + tot_exptime * u.s
                
                target_tmp['exec_start'] = obs_start_time.datetime
                target_tmp['exec_end'] = end_time.datetime
                if score is not None:
                    target_tmp['score'] = '%.2f'%score
                else:
                    target_tmp['score'] = score
                target_tmp['scheduled'] = True
                if isinstance(target_tmp, pd.core.series.Series):
                    target_df = target_tmp.to_frame().transpose()
                    target_df['scheduled'] = target_df['scheduled'].astype(bool, copy = False)
                else:
                    target_df = target_tmp
                output_table = pd.concat([output_table, target_df], axis = 0)
                return output_table, end_time
            
            def update_status_obs_tbl(obs_tbl,
                                    target):
                updated_tbl = obs_tbl.copy()
                idx = (updated_tbl['id'] == target['id'])
                updated_tbl.loc[idx, 'scheduled'] = True
                return updated_tbl
            
            def transitioner(name : str,
                            exptime : float or str,
                            count : float or str,
                            filter_ : str = '',
                            id_ : str = None
                            ):
                
                if id_ == None:
                    id_ = uuid.uuid4().hex
                transitioner = pd.DataFrame(data = [[str(name), str(exptime), str(count), str(filter_), str(id_), bool(True)]], columns = ['obj', 'exptime', 'count', 'filter', 'id', 'scheduled'])
                return transitioner.iloc[0]

            if ((len(schedule.unscheduled) > 0) | (len(scheduled_tbl) >0)):
                # Observation
                time_observation = schedule.obsnight.obs_start + delay_minute * u.minute
                autofocus_time = time_observation
                if autofocus_at_init:
                    autofocus_time = time_observation - 5 * u.hour
                telescope_altaz = None                    
                target_best = schedule.observable.iloc[0] 
                
                for i in range(n_target_for_each_timeslot):
                    time_observation = schedule.obsnight.obs_start  
                    while (time_observation < schedule.obsnight.obs_end):
                        # Autofocus
                        if autofocus:
                            remain_targets = schedule.observable[schedule.observable['scheduled'] == False]
                            if (Time(time_observation) - autofocus_time > autofocus_interval_in_hour * u.hour) & (duplicate_when_empty):
                                target_autofocus = transitioner(name = 'autofocus', exptime = self.config['OVERHEAD_TIME_AUTOFOCUS'], count = 1)
                                schedule.scheduled, time_observation = insert_schedule(obs_start_time = time_observation, scheduled_table=schedule.scheduled, target = target_autofocus, score = 1)
                                autofocus_time = time_observation
                                
                            elif (Time(time_observation) - autofocus_time > autofocus_interval_in_hour * u.hour) & (not duplicate_when_empty) & (len(remain_targets) > 0):
                                target_autofocus = transitioner(name = 'autofocus', exptime = self.config['OVERHEAD_TIME_AUTOFOCUS'], count = 1)
                                schedule.scheduled, time_observation = insert_schedule(obs_start_time = time_observation, scheduled_table=schedule.scheduled, target = target_autofocus, score = 1)
                                autofocus_time = time_observation
                                
                        immediate_targets = Table()
                        if len(scheduled_tbl) > 0:
                            immediate_targets = scheduled_tbl[(np.abs(scheduled_tbl['obstime_start'] - Time(time_observation).datetime) < datetime.timedelta(hours =2))]
                            immediate_targets = immediate_targets.sort_values(by ='obstime_start', ascending = True)
                        if len(immediate_targets) > 0:
                            immediate_target = immediate_targets.iloc[0]
                            #time_observation = immediate_target['obstime_start']
                            all_score = self.scorer(obs_tbl = immediate_targets, utctime = Time(time_observation), telescope_altaz = telescope_altaz, duplicate = False)
                            if np.sum(all_score) > 0:
                                target_best = immediate_target
                                time_observation = target_best['obstime_start']
                                target_best_score = all_score.iloc[0]
                                target = mainTarget(name_telescope = self.name_telescope, name_project = self.name_project, observer = self.observer, target_ra = target_best['ra'], target_dec = target_best['dec'])
                                telescope_altaz = target.altaz(utctimes = time_observation)
                                schedule.scheduled, time_observation = insert_schedule(obs_start_time = Time(time_observation), scheduled_table=schedule.scheduled, target = target_best, score = target_best_score)
                                target_overhead = transitioner(name = f'overhead_{target_best["obj"]}', exptime = calculate_tot_overheadtime(target_best), count = 1)
                                schedule.scheduled, time_observation = insert_schedule(obs_start_time = Time(time_observation), scheduled_table=schedule.scheduled, target = target_overhead)
                                schedule.observable = update_status_obs_tbl(schedule.observable, target = target_best)
                                schedule.all = update_status_obs_tbl(schedule.all, target = target_best)
                                scheduled_tbl = scheduled_tbl[~(scheduled_tbl['id'] == target_best['id'])]
                            else:
                                target_best = transitioner(name = f'empty_target', exptime = 600, count = 1)
                                schedule.scheduled, time_observation = insert_schedule(obs_start_time = Time(time_observation), scheduled_table=schedule.scheduled, target = target_best)
                        else:
                            # Target selection
                            all_score = self.scorer(obs_tbl = schedule.observable, utctime = Time(time_observation), telescope_altaz = telescope_altaz, duplicate = False)
                            if np.sum(all_score) > 0:
                                best_idx = np.argmax(all_score)
                                target_best = schedule.observable.iloc[best_idx] 
                                target_best_score = all_score.iloc[best_idx]
                                target = mainTarget(name_telescope = self.name_telescope, name_project = self.name_project, observer = self.observer, target_ra = target_best['ra'], target_dec = target_best['dec'])
                                telescope_altaz = target.altaz(utctimes = time_observation)
                                schedule.scheduled, time_observation = insert_schedule(obs_start_time = Time(time_observation), scheduled_table=schedule.scheduled, target = target_best, score = target_best_score)
                                target_overhead = transitioner(name = f'overhead_{target_best["obj"]}', exptime = calculate_tot_overheadtime(target_best), count = 1)
                                schedule.scheduled, time_observation = insert_schedule(obs_start_time = Time(time_observation), scheduled_table=schedule.scheduled, target = target_overhead)
                                schedule.observable = update_status_obs_tbl(schedule.observable, target = target_best)
                                schedule.all = update_status_obs_tbl(schedule.all, target = target_best)
                            else:
                                if duplicate_when_empty:
                                    all_targets = schedule.scheduled.dropna(subset = ['score'])
                                    last_targets = all_targets[-3:]
                                    duplicate_target = schedule.observable[~schedule.observable['id'].isin(last_targets['id'])]
                                    #previous_idx = (schedule.observable['id'] == target_best['id'])
                                    #extended_target = schedule.observable[~previous_idx]
                                    all_score = self.scorer(obs_tbl = duplicate_target, utctime = Time(time_observation), telescope_altaz = telescope_altaz, duplicate = duplicate_when_empty)
                                    if np.sum(all_score) > 0:
                                        best_idx = np.argmax(all_score)
                                        target_best = duplicate_target.iloc[best_idx] 
                                        target_best_score = all_score.iloc[best_idx]
                                        schedule.scheduled, time_observation = insert_schedule(obs_start_time = Time(time_observation), scheduled_table=schedule.scheduled, target = target_best, score = target_best_score)

                                    else:
                                        target_best = transitioner(name = f'empty_target', exptime = 600, count = 1)
                                        schedule.scheduled, time_observation = insert_schedule(obs_start_time = Time(time_observation), scheduled_table=schedule.scheduled, target = target_best)
                                else:
                                    target_best = transitioner(name = f'empty_target', exptime = 600, count = 1)
                                    schedule.scheduled, time_observation = insert_schedule(obs_start_time = Time(time_observation), scheduled_table=schedule.scheduled, target = target_best)

            
                if len(schedule.scheduled) > 0:    
                    schedule.scheduled = schedule.scheduled.dropna(subset = ['score'])
            else:
                pass
            
            schedule.all = Table().from_pandas(schedule.all)    
            schedule.observable = Table().from_pandas(schedule.observable)
            schedule.scheduled = Table().from_pandas(schedule.scheduled)
            schedule.scheduled.remove_column('coord.ra')
            schedule.scheduled.remove_column('coord.dec')
        else:
            schedule.all = Table().from_pandas(schedule.all)    
            schedule.observable = Table()
            schedule.scheduled = schedule.observable.copy()[0:0]
        return schedule
        
    def scorer(self,
                obs_tbl : Table,
                utctime : datetime or Time,
                telescope_altaz : SkyCoord = None,
                duplicate : bool = False
                ):
        utctime = Time(utctime)
        all_ra = np.array(obs_tbl['ra'])
        all_dec = np.array(obs_tbl['dec'])
        all_mainTarget = mainTarget(name_telescope = self.name_telescope, name_project =self.name_project, observer = self.observer, target_ra = all_ra, target_dec = all_dec)

        all_target_altaz = all_mainTarget.altaz(utctimes = utctime)
        all_target_alt = all_target_altaz.alt.value
        all_target_weight = obs_tbl['weight'].astype(float)
        all_target_hourangle_tmp = utctime.datetime - obs_tbl['transit']
        all_target_hourangle_sec = np.array([target_hourangle.total_seconds() for target_hourangle in all_target_hourangle_tmp])
        if not telescope_altaz == None:
            target_sepdist = SkyCoord.separation(telescope_altaz, all_target_altaz).value
        
        score = np.ones(len(obs_tbl))
        # constraints
        constraint_moonsep = obs_tbl['moonsep'] > self.config['TARGET_MOONSEP']
        score *= constraint_moonsep
        
        constraint_altitude = all_target_alt > self.config['TARGET_MINALT']
        score *= constraint_altitude
        
        constraint_altitude = all_target_alt < self.config['TARGET_MAXALT']
        score *= constraint_altitude
        
        constraint_meridian = (np.abs(all_target_hourangle_sec) > self.config['TARGET_MERIDIAN_SEC']) #(-self.config['TARGET_MERIDIAN_SEC'] > all_target_hourangle) |  ( all_target_hourangle > 0)
        #score *= constraint_meridian
        
        if not telescope_altaz == None:
            constraint_maxslew = target_sepdist < self.config['TARGET_MAXSLEW']
            score *= constraint_maxslew

        if not duplicate:
            constraint_duplicate = np.array(~obs_tbl['scheduled'])
            score *= constraint_duplicate    
        
        # score
        all_target_alt = np.array([0 if target_alt <= 0 else target_alt for target_alt in all_target_alt])
        score_relative_alt = self.config['TARGET_WEIGHT_RELATIVE_ALT'] * (all_target_alt) / (np.abs(obs_tbl['maxalt']) + 0.1)
        score_absolute_alt = self.config['TARGET_WEIGHT_ABSOLUTE_ALT'] * (all_target_alt) / (90)
        
        max_weight = np.max(all_target_weight)
        score_weight = self.config['TARGET_WEIGHT_PRIORITY']* (all_target_weight / max_weight)
        if self.name_project == 'IMSNG':
            score_weight = self.config['TARGET_WEIGHT_PRIORITY']* (1- (all_target_weight / max_weight))
        
        if not telescope_altaz == None:
            score_slewtime = self.config['TARGET_WEIGHT_SLEW'] * (target_sepdist / self.config['OVERHEAD_SPEED_SLEWING']) / np.max((target_sepdist / self.config['OVERHEAD_SPEED_SLEWING']))
            score_all = (score_relative_alt + score_absolute_alt + score_weight + score_slewtime) / (self.config['TARGET_WEIGHT_RELATIVE_ALT'] + self.config['TARGET_WEIGHT_ABSOLUTE_ALT'] + self.config['TARGET_WEIGHT_PRIORITY'] + self.config['TARGET_WEIGHT_SLEW']) 
        else:
            score_all = (score_relative_alt + score_absolute_alt + score_weight) / (self.config['TARGET_WEIGHT_RELATIVE_ALT'] + self.config['TARGET_WEIGHT_ABSOLUTE_ALT'] + self.config['TARGET_WEIGHT_PRIORITY']) 
        score *= score_all
        return score

    def show(self,
                scheduled_table : Table,
                save : bool = False,
                filename_prefix : str = 'ToO_',
                savepath = './ACPscript/'
                ):
        
        now = Time.now()
        
        if len(scheduled_table) > 0:
            scheduled_table_pd = scheduled_table.to_pandas()
            target_table_pd = scheduled_table_pd.dropna(subset = ['ra'])
            calib_table_pd = scheduled_table_pd[scheduled_table_pd['ra'].isna()]
            target_table = Table().from_pandas(target_table_pd)
            calib_table = Table().from_pandas(calib_table_pd)
            n_obj = len(set(target_table['obj']))
            n_obs = len(set(target_table['exec_start']))
            
            time_astronomical_twilight = self.obsnight.sunset_astro
            time_astronomical_dawn = self.obsnight.sunrise_astro
            time_civil_twilight = self.obsnight.sunset_prepare
            time_civil_dawn = self.obsnight.sunrise_prepare
            time_obs_start = self.obsnight.obs_start
            time_obs_end = self.obsnight.obs_end

            time_range_start, time_range_end = time_civil_twilight - 2*u.hour, time_civil_dawn + 2*u.hour
            time_axis = np.arange(time_range_start.datetime, time_range_end.datetime, datetime.timedelta(minutes = 5))
            moon_altaz = self.observer.moon_altaz(time_axis)
            sun_altaz = self.observer.sun_altaz(time_axis)
            
            plt.figure(dpi = 300, figsize = (10, 4))
            plt.title(f'Observing plan[{self.name_telescope}] on {time_obs_end.datetime.strftime("%Y-%m-%d")}')
            if len(target_table) > 0:
                for target in target_table:
                    target_tcspy = mainTarget(name_telescope = self.name_telescope, name_project = self.name_project, observer = self.observer, target_ra = target['ra'], target_dec = target['dec'])
                    target_altaz = target_tcspy.altaz(time_axis)
                    obs_time_axis = np.arange(target['exec_start'].datetime, target['exec_end'].datetime, datetime.timedelta(minutes = 1))
                    obs_target_altaz = target_tcspy.altaz(obs_time_axis)
                    plt.plot(target_altaz.obstime.datetime, target_altaz.alt.value, c='k', linewidth = 0.3)
                    plt.plot(obs_target_altaz.obstime.datetime, obs_target_altaz.alt.value, c='r')
            if (now.datetime < time_range_end + datetime.timedelta(hours = 3)) & (now.datetime > time_range_start - datetime.timedelta(hours = 3)) * (~self._entire_night):
                plt.axvline(now.datetime, linestyle = '--', c='r', label = 'Now')
            plt.plot(moon_altaz.obstime.datetime, moon_altaz.alt.value, c = 'b', label ='Moon', linestyle = ':')
            plt.plot(sun_altaz.obstime.datetime, sun_altaz.alt.value, c = 'r', label = 'Sun', linestyle = ':')

            plt.fill_betweenx([10,90], time_astronomical_twilight.datetime, time_astronomical_dawn.datetime, alpha = 0.1, color='k')
            plt.fill_betweenx([10,90], time_civil_twilight.datetime, time_civil_dawn.datetime, alpha = 0.05, color='k')
            plt.axvline(x=time_astronomical_twilight.datetime, linestyle = '-', c='k', linewidth = 0.5)
            plt.axvline(x=time_astronomical_dawn.datetime, linestyle = '-', c='k', linewidth = 0.5)
            plt.axvline(x=time_civil_dawn.datetime, linestyle = '--', c='k', linewidth = 0.5)
            plt.axvline(x=time_civil_twilight.datetime, linestyle = '--', c='k', linewidth = 0.5)
            plt.axhline(y=self.config['TARGET_MINALT'], linestyle = '--', c='r', linewidth = 1)
            plt.text(time_civil_twilight.datetime-datetime.timedelta(minutes=140), 85, f'N_observation[N_target] = {(n_obs)}[{n_obj}]', fontsize = 10)
            plt.text(time_astronomical_twilight.datetime, 92, 'Ast.S.set', fontsize = 10)
            plt.text(time_civil_twilight.datetime-datetime.timedelta(minutes=00), 92, 'S.set', fontsize = 10)
            plt.text(time_civil_dawn.datetime-datetime.timedelta(minutes=00), 92, 'S.rise', fontsize = 10)
            #plt.text(time_civil_dawn.datetime-datetime.timedelta(minutes=00), 40, f'Input time\n {now.datetime.strftime("%m-%d-%H")}', fontsize = 10)
            plt.xlim(time_range_start.datetime - datetime.timedelta(hours = 1), time_range_end.datetime + datetime.timedelta(hours = 1))
            if len(calib_table) > 0:
                for target in calib_table:
                    plt.fill_betweenx([10,90], target['exec_start'].datetime, target['exec_end'].datetime, alpha = 0.1, color ='r')
                plt.fill_betweenx([-30,-31], target['exec_start'].datetime, target['exec_end'].datetime, alpha = 0.1, color ='r', label = 'Calib, Focus')
            plt.ylim(10, 90)
            plt.legend(loc = 1)
            plt.xlabel('UTC [mm-dd hh]')
            plt.ylabel('Altitude [deg]')
            plt.grid()
            plt.xticks(rotation = 45)
            if save:
                os.makedirs(savepath, exist_ok = True)
                filename = f'{savepath}{filename_prefix}ObservingSchedule_{(self.obsnight.midnight.datetime).strftime("%y%m%d")}_{self.name_telescope}.png'
                plt.savefig(filename, bbox_inches = 'tight')
                print(f'[show] SAVED! filename : {filename}')
        else:
            print('[show] FAILED! No observable Target.')    
            
    def _set_target(self):
        # Set all target & observable tatget objects
        # Update maximum altitude and transit time of the targets (55s/24000targets)
        
        class alltarget: pass
        target_db = self._target_db.copy()
        if isinstance(target_db, pd.core.frame.DataFrame):
            target_table = Table().from_pandas(target_db)
        elif isinstance(target_db, Table):
            target_table= target_db.copy()
        elif isinstance(target_db, dict):
            header = target_db['header']
            values = target_db['value']
            data_pd = pd.DataFrame(values, columns = header)
            target_table = Table().from_pandas(data_pd)
        else:
            raise ValueError("variable 'target_db' must be pd.Dataframe, astropy.table.Table, dict format")
        
        default_values = dict(binning = 1, 
                              weight = 1)
        
        for column in ['obj', 'ra', 'dec', 'filter', 'exptime', 'count']:
            if not column in target_table.keys():
                raise ValueError(f'"{column}" does not exist in the data')
        for column in ['binning']:
            if not column in target_table.keys():
                target_table[column] = default_values[column]
        for column in ['weight']:
            if not column in target_table.keys():
                try:
                    target_table[column] = target_table['priority']
                except:
                    target_table[column] = default_values[column]
                
        for column in ['exptime', 'count', 'binning']:
            col_valuelist = []
            for target in target_table:
                len_filt = len(str(target['filter']).split(','))
                len_col = len(str(target[column]).split(','))
                if len_filt != len_col:
                    if len_filt == 0:
                        value = str(default_values[column]).split(',')[0]
                    else:
                        value = str(target[column]).split(',')[0]
                    col_value = ','.join([value] * len_filt)
                    col_valuelist.append(col_value)
                else:
                    col_value = str(target[column])
                    col_valuelist.append(col_value)
            target_table[column] = col_valuelist
                
        try:        
            coord = self._match_coord_format(target_table['ra'], target_table['dec'])            
        except:
            coord = self._match_coord_format(target_table['ra_hms'], target_table['dec_dms'])   
        ra_string, dec_string = np.array(self._match_coord_to_string(coord))
        target_table['coord'] = coord
        target_table['ra'] = ['%.5f'%ra for ra in target_table['coord'].ra.value]
        target_table['dec'] = ['%.5f'%dec for dec in target_table['coord'].dec.value]
        target_table['ra_hms'] = ra_string
        target_table['dec_dms'] = dec_string
        
        transittime, maxalt = self._get_maxaltaz(target_table)
        moonsep = self._get_moonsep(target_table)
        target_table['maxalt'] = maxalt
        target_table['transit'] = transittime
        target_table['moonsep'] = moonsep
        if not 'scheduled' in target_table.colnames:
            target_table['scheduled'] = False     
        else:
            try:
                def str_to_bool(value):
                    return str(value).lower() == 'true'
                scheduled_array = np.array([str_to_bool(value) for value in target_table['scheduled']])
                target_table.remove_column('scheduled')
                target_table['scheduled'] = scheduled_array
            except:
                target_table['scheduled'] = target_table['scheduled'].astype(bool)
            
        target_table['id'] = [uuid.uuid4().hex for i in range(len(target_table))]
        alltarget.observable = self._get_target_observable(obs_tbl = target_table, fraction_observable = 0.1)
        alltarget.all = target_table
        return alltarget
    
    def _get_moonsep(self,
                        obs_tbl : Table):
        all_coords = obs_tbl['coord']
        moon_coord = SkyCoord(ra =self.obsinfo.moon_radec.ra.value, dec = self.obsinfo.moon_radec.dec.value, unit = 'deg')
        moonsep = np.array(SkyCoord.separation(all_coords, moon_coord).value).round(2)
        return moonsep
    
    def _get_maxaltaz(self,
                        obs_tbl : Table):
        all_ra = np.array(obs_tbl['ra'])
        all_dec = np.array(obs_tbl['dec'])
        all_coords = obs_tbl['coord']
        all_mainTarget = mainTarget(name_telescope = self.name_telescope, name_project = self.name_project, observer = self.observer, target_ra = all_ra, target_dec = all_dec)
        all_time_hourangle = all_mainTarget.hourangle(self.obsnight.midnight)
        all_hourangle_converted = [hourangle if (hourangle -12 < 0) else hourangle-24 for hourangle in all_time_hourangle.value]
        all_target_altaz_at_sunset = all_mainTarget.altaz(utctimes=self.obsnight.sunset_astro)
        all_target_altaz_at_sunrise = all_mainTarget.altaz(utctimes=self.obsnight.sunrise_astro)
        all_transittime = self.obsnight.midnight - all_hourangle_converted * u.hour
        all_maxalt = []
        for i, target_info in enumerate(zip(all_transittime, all_coords)):
            target_time_transit, target_coord = target_info
            if (target_time_transit > self.obsnight.sunset_astro) & (target_time_transit < self.obsnight.sunrise_astro):
                maxaltaz = self.obsinfo.observer_astroplan.altaz(target_time_transit, target = target_coord)
                maxalt = np.round(maxaltaz.alt.value,2)
            else:
                sunset_alt = all_target_altaz_at_sunset[i].alt.value
                sunrise_alt = all_target_altaz_at_sunrise[i].alt.value
                maxalt = np.round(np.max([sunset_alt, sunrise_alt]),2)
            all_maxalt.append(maxalt)
        return all_transittime, all_maxalt
        
    def _get_target_observable(self,
                                obs_tbl : Table,
                                fraction_observable : float = 0.1):
        # obs_tbl must have ra, dec, obj column
        all_targets = obs_tbl['coord']
        constraints = self.constraint.astroplan
        observer = self.obsinfo.observer_astroplan
        night_start = self.obsnight.sunset_astro
        night_end = self.obsnight.sunrise_astro
        observability_tbl = observability_table(constraints = constraints, observer = observer, targets = all_targets , time_range = [night_start, night_end], time_grid_resolution = 20 * u.minute)
        obs_tbl['fraction_obs'] = ['%.2f'%fraction for fraction in observability_tbl['fraction of time observable']]
        key = observability_tbl['fraction of time observable'] > fraction_observable
        obs_tbl = obs_tbl[key]
        
        return obs_tbl
        
    def _set_obs_info(self):
        class info: pass
        time_input = Time(self._date)
        info.moon_phase = self.observer.moon_phase(time_input)
        info.moon_radec = self.observer.moon_radec(time_input)
        info.sun_radec = self.observer.sun_radec(time_input)
        info.observer_info = self.observer.get_info()
        info.observer_astroplan = self.observer._observer
        info.is_night = self.observer.is_night(time_input)
        return info
        
    def _set_obs_night(self,
                        horizon_prepare : float = -5,
                        horizon_astro : float = -15):
        class night: pass
        time_input = Time(self._date)
        night.sunrise_prepare = self.observer.tonight(time = time_input, horizon = horizon_prepare)[1]
        night.sunset_prepare = self.observer.sun_settime(night.sunrise_prepare, mode = 'previous', horizon= horizon_prepare)
        night.sunrise_astro = self.observer.sun_risetime(night.sunrise_prepare, mode = 'previous', horizon= horizon_astro)
        night.sunset_astro = self.observer.sun_settime(night.sunrise_prepare, mode = 'previous', horizon= horizon_astro)
        night.midnight = Time((night.sunset_astro.jd + night.sunrise_astro.jd)/2, format = 'jd')
        if self._entire_night:
            night.obs_start = night.sunset_astro
        else:
            if (time_input < night.sunrise_astro) & (time_input > night.sunset_astro):
                night.obs_start = time_input
            else:
                night.obs_start = night.sunset_astro
        night.obs_end = night.sunrise_astro
        
        return night
    
    def _set_constrints(self):
        class constraint: pass
        constraint_astroplan = []
        if (self.config['TARGET_MINALT'] != None) & (self.config['TARGET_MAXALT'] != None):
            constraint_altitude = AltitudeConstraint(min = self.config['TARGET_MINALT'] * u.deg, max = self.config['TARGET_MAXALT'] * u.deg, boolean_constraint = False)
            constraint_astroplan.append(constraint_altitude)
            constraint.maxalt = self.config['TARGET_MINALT']
            constraint.minalt = self.config['TARGET_MAXALT']
        if self.config['TARGET_MAX_SUNALT'] != None:
            constraint_atnight = AtNightConstraint(max_solar_altitude = self.config['TARGET_MAX_SUNALT'] * u.deg)
            constraint_astroplan.append(constraint_atnight)
            constraint.sun_maxalt = self.config['TARGET_MAX_SUNALT']
        if self.config['TARGET_MOONSEP'] != None:
            constraint_gallatitude = MoonSeparationConstraint(min = self.config['TARGET_MOONSEP'] * u.deg, max = None)
            constraint_astroplan.append(constraint_gallatitude)
            constraint.moon_separation = self.config['TARGET_MOONSEP']
        if self.config['TARGET_MAXAIRMASS'] != None:
            constraint_gallatitude = AirmassConstraint(min = 1, max = self.config['TARGET_MAXAIRMASS'], boolean_constraint = False)
            constraint_astroplan.append(constraint_gallatitude)
            constraint.airmass = self.config['TARGET_MAXAIRMASS']
        constraint.astroplan = constraint_astroplan
        return constraint

    def _match_coord_format(self, ra, dec):
        """
        Create a SkyCoord object from input coordinates ra and dec.

        Parameters:
        - ra (float, int, str): Right ascension in various formats.
        - dec (float, int, str): Declination in various formats.

        Returns:
        - skycoord (SkyCoord): SkyCoord object.

        Raises:
        - ValueError: If input format is not supported.
        """
        if isinstance(ra, (float, int, str)):
            ra = str(ra)
            dec = str(dec)
            if (':' in ra) and (':' in dec):
                skycoord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
            elif ('h' in ra) and ('d' in dec):
                skycoord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
            elif (' ' in ra) and (' ' in dec):
                skycoord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
            else:
                skycoord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
        else:
            try:
                ra0 = str(ra[0])
                dec0 = str(dec[0])
                if (':' in ra0) and (':' in dec0):
                    skycoord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
                elif ('h' in ra0) and ('d' in dec0):
                    skycoord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
                elif (' ' in ra0) and (' ' in dec0):
                    skycoord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
                else:
                    skycoord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
            except Exception as e:
                raise ValueError(f'Input format is not supported: {str(e)}')
        return skycoord
        
    def _match_coord_to_string(self, 
                                coords : SkyCoord):
        # This is more than 50x faster than coords.to_string(sep = ':')
        def string_ra(ra):
            ra_hour = ra.hms.h
            ra_minute = ra.hms.m
            ra_seconds = ra.hms.s
            ra_string = '%02d:%02d:%02d'%(ra_hour, ra_minute, ra_seconds)
            return ra_string
            
        def string_dec(dec):
            dec_deg = np.abs(dec.dms.d)
            dec_minute = np.abs(dec.dms.m)
            dec_seconds = np.abs(dec.dms.s)
            if dec.value < 0:
                dec_string = '-%02d:%02d:%02d'%(dec_deg, dec_minute, dec_seconds)
            else:
                dec_string = '%02d:%02d:%02d'%(dec_deg, dec_minute, dec_seconds)
            return dec_string
        ra_string = list(map(string_ra, coords.ra))
        dec_string = list(map(string_dec, coords.dec))
        return ra_string, dec_string


class ScriptMaker(mainConfig):
    """
    Support module to make a script based on the ObsScheduler object information
    - History
    (23.05.02) ""Version 3.0"", Running time decreased.
        
    ==========
    Parameters
    ==========
    1. ObsScheduler : ObsScheduler
        ObsScheduler instance to make a script. 
        If No schedule is defined (ObsScheduler.schedule == None), 
        all the methods in ScriptMaker will make a new schedule.
        
    =======
    Methods
    =======
    1. write_rts(n_target_for_each_timeslot : int, **kwargs) -> None
        Write the RTS script with the scheduled table(ObsScheduler.schedule.scheduled)
        - n_target_for_each_timeslot : int = if 2, 2 times of the scheduled target will be written 
        
    2. write_ACPscript_RASA36(**kwargs) -> None
        Write the ACP script for RASA36 telescope
        
    3. write_ACPscript_LSGT(period_script: int, **kwargs)
        Write the ACP script for LSGT telescope
        - period_script : int = the script will be divided by period_script (hour)
    
    4. write_ACPscript_KCT(**kwargs) -> None
        Write the ACP script for KCT telescope
    
    5. show(**kwargs) -> None
        Visualize the observing plan with the scheduled table 
    
    6. write_log(**kwargs) -> log_table : Table
        Write log file for the scheduled table. 
        All scheduled targets must be included in the log file.
        
        - n_target : int = the length of the target in the log file
        - sort_keyword : str = the keyword to sort the log 
        - return_ : bool = if True, return log_table
    
    ==========
    Properties
    ==========
    1. result : ObsScheduler object (See properties of ObsScheduler class)
    """
    
    def __init__(self, 
                 ObsScheduler):
        self.obsdata = ObsScheduler
    
    @property 
    def result(self):
        return self.obsdata.schedule
    
    def write_rts(self,
                  filename_prefix : str = '',
                  savepath = './rts/',
                  n_target_for_each_timeslot : int = 2,
                  **kwargs):
        
        def dt_to_form_hoursec(dt):
            form ='{:02}:{:02}'.format(dt.hour, dt.minute)
            return form
        
        def ut_to_lt(utctime):
            return self.obsdata.observer.localtime(utctime)
        
        os.makedirs(savepath, exist_ok = True)
        
        if self.obsdata.schedule == None:
            self.obsdata.schedule = self.obsdata.scheduler(obs_tbl = self.obsdata.target.observable, n_target_for_each_timeslot = n_target_for_each_timeslot)#, **kwargs)
        #else:
        #    self.obsdata.schedule = self.obsdata.scheduler(obs_tbl = self.obsdata.schedule.observable, n_target_for_each_timeslot = n_target_for_each_timeslot)#, **kwargs)
        
        if len(self.obsdata.schedule.scheduled) > 0:
            scheduled_tbl_pd = self.obsdata.schedule.scheduled.to_pandas()
            scheduled_tbl = scheduled_tbl_pd.dropna(subset = ['ra']).copy()
            # Target (rts)
            obj_target = np.array([mainTarget(name_telescope= self.obsdata.name_telescope, name_project = self.obsdata.name_project, observer = self.obsdata.observer, target_ra = target['ra'], target_dec = target['dec']) for idx, target in scheduled_tbl.iterrows()])
            scheduled_tbl['mainTarget'] = obj_target
            
            risetime = []
            settime = []
            transittime = []
            for i, target in scheduled_tbl.iterrows():
                transittime.append(dt_to_form_hoursec(ut_to_lt(target['transit'])))
                try:
                    risetime.append(dt_to_form_hoursec(ut_to_lt(target['mainTarget'].risetime(target['transit'], n_grid_points = 100, mode = 'previous').datetime)))
                except:
                    risetime.append(dt_to_form_hoursec(ut_to_lt(self.obsdata.obsnight.sunset_astro.datetime)))
                try:
                    settime.append(dt_to_form_hoursec(ut_to_lt(target['mainTarget'].settime(target['transit'], n_grid_points = 100, mode = 'next').datetime)))
                except:
                    settime.append(dt_to_form_hoursec(ut_to_lt(self.obsdata.obsnight.sunrise_astro.datetime)))
            scheduled_tbl['transit(LT)'] = transittime
            scheduled_tbl['rise(LT)'] = risetime
            scheduled_tbl['set(LT)'] = settime
            # Moon
            moon_radec = self.obsdata.obsinfo.moon_radec
            moonphase = round(self.obsdata.obsinfo.moon_phase,2)
            moon_ra_str, moon_dec_str = moon_radec.to_string('hmsdms', sep = ':', precision = 0).split(' ')
            moon_set_sepration = int(self.obsdata.config['TARGET_MOONSEP'])
            # Sun   
            sunrisetime = dt_to_form_hoursec(ut_to_lt(self.obsdata.obsnight.sunrise_astro.datetime))
            sunsettime = dt_to_form_hoursec(ut_to_lt(self.obsdata.obsnight.sunset_astro.datetime))
            tonight = ut_to_lt(self.obsdata.obsnight.midnight.datetime).strftime('%Y-%m-%d')
            filename = f'{savepath}{filename_prefix}rts_{(ut_to_lt(self.obsdata.obsnight.midnight.datetime)).strftime("%y%m%d")}_{self.obsdata.name_telescope}.txt'
            # Target table
            rts_table = Table()
            rts_table['name'] = scheduled_tbl['obj']
            rts_table['ra'] = scheduled_tbl['ra_hms']
            rts_table['dec'] = scheduled_tbl['dec_dms']
            rts_table['rise(LT)'] = scheduled_tbl['rise(LT)']
            rts_table['transit(LT)'] = scheduled_tbl['transit(LT)']
            rts_table['set(LT)'] = scheduled_tbl['set(LT)']
            rts_table['moon_dist(deg)'] = scheduled_tbl['moonsep']
            rts_table['weight'] = scheduled_tbl['weight']
            rts_table['note'] = scheduled_tbl['note']
            #rts_table['filter'] = obs_info['filter']
            #rts_table['exptime'] = obs_info['exptime']
            #rts_table['count'] = obs_info['count']
            #rts_table['binning'] = obs_info['binning']
            # Sorting
            rts_table['ra_sort'] = scheduled_tbl['ra']
            rts_table.sort('ra_sort')
            rts_table.remove_column('ra_sort')
            rts_table.write(filename, format = 'ascii.fixed_width', overwrite = True)
            print(f'[write_rts] SAVED! filename : {filename}')
            
            # add info
            f = open(filename, 'r')
            original_data = f.readlines()
            with open(filename, "w") as f:
                f.write('#\tObserving Date = '+tonight+'\n')
                f.write('#\tObservatory = '+self.obsdata.name_telescope+'\n')
                f.write('#\t-15 deg sunset = '+str(sunsettime)+'\n')
                f.write('#\t-15 deg sunrise = '+str(sunrisetime)+'\n')
                f.write('#\tMoon ra, dec = '+ str(moon_ra_str) +'   '+ str(moon_dec_str)+'\n')
                f.write('#\tMoon phase = '+str(moonphase)+'\n')
                f.write('#\tMoon seperation limit = '+str(moon_set_sepration)+'\n')
                f.write('#\tAltitude limit = '+str(self.obsdata.config['TARGET_MINALT'])+'\n')
                f.write('='*110+'\n')
                for line in original_data:
                    f.write(line)
                f.write('='*110+'\n')
        else:
            print(f'[write_rts] FAILED! No observable target!')
            
    def write_ACPscript_RASA36(self,
                               duplicate_when_empty : bool = False,
                               n_target_for_each_timeslot : int = 1,
                               delay_minute : int = 10,
                               
                               ccdcool : bool = True,
                               ccdtemp : float = -10,
                               shutdown : bool = False,
                               
                               autofocus : bool = True,
                               autofocus_at_init : bool = False,
                               autofocus_interval_in_hour : float = 3,
                               
                               bias : bool = False,
                               bias_count : int = 9,
                               dark : bool = False,
                               dark_count : int = 9,
                               dark_exptime : float = 60,
                               flat_dusk : bool = False,
                               flat_dawn : bool = False,
                               
                               filename_prefix : str = 'ToO_',
                               savepath = './ACPscript/'
                               ):
        
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form
        
        
        os.makedirs(savepath, exist_ok = True)
        if self.obsdata.schedule == None:
            self.obsdata.schedule = self.obsdata.scheduler(obs_tbl = self.obsdata.target.observable, n_target_for_each_timeslot = n_target_for_each_timeslot, autofocus = autofocus, autofocus_at_init = autofocus_at_init, autofocus_interval_in_hour = autofocus_interval_in_hour, duplicate_when_empty = duplicate_when_empty, delay_minute = delay_minute)

        scheduled_table = self.obsdata.schedule.scheduled
        if len(scheduled_table) > 0:
            
            # Set observing time 
            time_civil_twilight = self.obsdata.obsnight.sunset_prepare
            time_civil_dawn = self.obsdata.obsnight.sunrise_prepare
            time_astronomical_twilight = self.obsdata.obsnight.sunset_astro
            time_astronomical_dawn = self.obsdata.obsnight.sunrise_astro
            time_obs_start = self.obsdata.obsnight.obs_start
            time_obs_end = self.obsdata.obsnight.obs_end
            
            time_civil_twilight_form = dt_to_form(time_civil_twilight.datetime)
            time_civil_dawn_form = dt_to_form(time_civil_dawn.datetime)
            time_astronomical_twilight_form = dt_to_form(time_astronomical_twilight.datetime)
            time_astronomical_dawn_form = dt_to_form(time_astronomical_dawn.datetime)
            time_obs_start_form = dt_to_form(time_obs_start.datetime)
            time_obs_end_form = dt_to_form(time_obs_end.datetime)
            filename = f'{savepath}{filename_prefix}ACP_{(self.obsdata.obsnight.midnight.datetime).strftime("%y%m%d_%H%M")}_{self.obsdata.name_telescope}.txt'
        

            with open(filename, 'w') as f:
                
                # Preparation
                f.write('; Prepare for the observation (RASA36)\n')
                f.write('#WAITUNTIL 1, {}\n'.format(time_civil_twilight_form))
                
                # Cooling
                if ccdcool:
                    f.write('\n; Cooling\n')
                    f.write('#CHILL {}, {}\n'.format(ccdtemp, 1.0))
                
                # Calibration frames
                f.write('\n; Calibration frames\n')
                if bias:
                    f.write("#COUNT {}\n#BIAS\n".format(bias_count)) # bias
                    
                if dark:
                    f.write("#COUNT {}\n#INTERVAL {}\n#DARK\n".format(dark_count, dark_exptime)) # dark
                
                if flat_dusk:
                    f.write('#DUSKFLATS\n') # duskflat
                

                # Image frames
                f.write('\n; Start of the evening\n')
                f.write('#WAITUNTIL 1, {}\n'.format(time_astronomical_twilight_form))
                
                f.write('\n; Targeting\n')
                target_table = Table().from_pandas(scheduled_table.to_pandas().dropna(subset = ['ra']))
                for target in target_table:
                    if target['obj'] == 'autofocus':
                        f.write('\n#AUTOFOCUS\n')
                    else:
                        name = target['obj']
                        ra = target['ra_hms']
                        dec = target['dec_dms']
                        exptime = target['exptime']
                        count = target['count']
                        binning = target['binning']
                        f.write('#COUNT {}\n'.format(count))
                        f.write('#INTERVAL {}\n'.format(exptime))
                        #f.write('#BINNING {}\n'.format(binning))
                        f.write('#NOPREVIEW \n')
                        f.write('#NOSOLVE\n')
                        #f.write('#NOPOINTING\n')
                        f.write('{}\t{}\t{}\n'.format(name,ra,dec))
                        f.write('\n')

                f.write('#QUITAT {}\n'.format(time_astronomical_dawn_form))
                
                if flat_dawn:
                    f.write('#DAWNFLATS\n') # dawnflat
                
                if shutdown:
                    f.write('\n; Closing\n')
                    f.write('#SHUTDOWN\n')
            print(f'[write_ACPscript_RASA36] SAVED! filename : {filename}')
        else:
            print('[write_ACPscript_RASA36] FAILED! No observable Target.')
        
    def write_ACPscript_LSGT(self,
                             period_script : int = 3, # hour
                             duplicate_when_empty= False,
                             n_target_for_each_timeslot : int = 1,
                             delay_minute : int = 10,
                             
                             ccdcool : bool = True,
                             ccdtemp : float = -10,
                             shutdown : bool = False,
                             
                             autofocus : bool = True,
                             autofocus_at_init : bool = False,
                             autofocus_interval_in_hour : float = 3,
                             
                             bias : bool = False,
                             bias_count : int = 9,
                             dark : bool = False,
                             dark_count : int = 9,
                             dark_exptime : float = 60,
                             flat_dusk : bool = False,
                             flat_dawn : bool = False,
                             
                             filename_prefix : str = 'ToO_',
                             savepath = './ACPscript/'
                             ):
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form
        
        def ut_to_lt(utctime):
            return self.obsdata.observer.localtime(utctime)
        
        os.makedirs(savepath, exist_ok = True)
        if self.obsdata.schedule == None:
            self.obsdata.schedule = self.obsdata.scheduler(obs_tbl = self.obsdata.target.observable, n_target_for_each_timeslot = n_target_for_each_timeslot, autofocus = autofocus, autofocus_at_init = autofocus_at_init, autofocus_interval_in_hour = autofocus_interval_in_hour, duplicate_when_empty = duplicate_when_empty, delay_minute = delay_minute)
        scheduled_table = self.obsdata.schedule.scheduled
        if len(scheduled_table) > 0:
            target_table = Table().from_pandas(scheduled_table.to_pandas())
            target_table = target_table[target_table['id'] != '']
            if len(target_table) > 0:
            
                # Set observing time 
                time_civil_twilight = self.obsdata.obsnight.sunset_prepare
                time_civil_dawn = self.obsdata.obsnight.sunrise_prepare
                time_astronomical_twilight = self.obsdata.obsnight.sunset_astro
                time_astronomical_dawn = self.obsdata.obsnight.sunrise_astro
                time_obs_start = self.obsdata.obsnight.obs_start
                time_obs_end = self.obsdata.obsnight.obs_end
                
                time_civil_twilight_form = dt_to_form(time_civil_twilight.datetime)
                time_civil_dawn_form = dt_to_form(time_civil_dawn.datetime)
                time_astronomical_twilight_form = dt_to_form(time_astronomical_twilight.datetime)
                time_astronomical_dawn_form = dt_to_form(time_astronomical_dawn.datetime)
                time_obs_start_form = dt_to_form(time_obs_start.datetime)
                time_obs_end_form = dt_to_form(time_obs_end.datetime)
                
                time_start = target_table[0]['exec_start']
                time_tot = (target_table[-1]['exec_end'] - target_table[0]['exec_start']).value * 24
                time_grid = Time(time_start) + period_script * np.arange(time_tot//period_script + 1) * u.hour
                
                for timecut in time_grid:
                    if timecut != time_grid[-1]:
                        idx = (timecut <= target_table['exec_start']) & ( target_table['exec_start'] < timecut + period_script * u.hour)
                        split_table = target_table[idx]
                        if len(split_table) > 0:
                            time_script = np.round((split_table[-1]['exec_end'] - split_table[0]['exec_start']).value * 24,1)
                            filename = f'{savepath}{filename_prefix}ACP_{(ut_to_lt(split_table[0]["obs_start"].datetime)).strftime("%y%m%d_%H%M")}_{time_script}Hour_{self.obsdata.name_telescope}.txt'
                            with open(filename, 'w') as f:
                                f.write('; Prepare for the observation (LSGT)\n')
                                f.write('\n; Targeting\n')
                                for target in split_table:
                                    name = target['obj']
                                    ra = target['ra_hms']
                                    dec = target['dec_dms']
                                    exptime = target['exptime']
                                    filter_ = target['filter']
                                    count = target['count']
                                    binning = target['binning']
                                    f.write('#skippreviews \n')
                                    f.write('#filteroffsets\n')
                                    f.write('#count {}\n'.format(count))
                                    f.write('#interval {}\n'.format(exptime))
                                    f.write('#binning {}\n'.format(binning))
                                    f.write('#filter {}\n'.format(filter_))
                                    f.write('{}\t{}\t{}\n'.format(name,ra,dec))
                                    f.write('\n')
                            print(f'[write_ACPscript_LSGT] SAVED! filename : {filename}')
                    else:
                        idx = (timecut < target_table['exec_start'])
                        split_table = target_table[idx]
                        if len(split_table) > 0:
                            time_script = np.round((split_table[-1]['exec_end'] - split_table[0]['exec_start']).value * 24,1)
                            filename = f'{savepath}{filename_prefix}ACP_{(ut_to_lt(split_table[0]["obs_start"].datetime)).strftime("%y%m%d_%H%M")}_{time_script}Hour_{self.obsdata.name_telescope}.txt'
                            with open(filename, 'w') as f:
                                f.write('; Prepare for the observation\n')
                                f.write('\n; Targeting\n')
                                for target in split_table:
                                    name = target['obj']
                                    ra = target['ra_hms']
                                    dec = target['dec_dms']
                                    exptime = target['exptime']
                                    filter_ = target['filter']
                                    count = target['count']
                                    binning = target['binning']
                                    f.write('#skippreviews \n')
                                    f.write('#filteroffsets\n')
                                    f.write('#count {}\n'.format(count))
                                    f.write('#interval {}\n'.format(exptime))
                                    f.write('#binning {}\n'.format(binning))
                                    f.write('#filter {}\n'.format(filter_))
                                    f.write('{}\t{}\t{}\n'.format(name,ra,dec))
                                    f.write('\n')
                                f.write('#QUITAT {}\n'.format(time_astronomical_dawn_form))
                            print(f'[write_ACPscript_LSGT] SAVED! filename : {filename}')
                    
            else:
                print('[write_ACPscript_LSGT] FAILED! No observable Target.')

    def write_ACPscript_KCT(self,
                            duplicate_when_empty : bool = False,
                            n_target_for_each_timeslot : int = 1,
                            delay_minute : int = 10,
                            
                            ccdcool : bool = True,
                            ccdtemp : float = -10,
                            shutdown : bool = True,
                            
                            autofocus : bool = True,
                            autofocus_at_init : bool = True,
                            autofocus_interval_in_hour : float = 3,
                            
                            bias : bool = False,
                            bias_count : int = 9,
                            dark : bool = False,
                            dark_count : int = 9,
                            dark_exptime : float = 60,
                            flat_dusk : bool = False,
                            flat_dawn : bool = False,
                            
                            filename_prefix : str = 'ToO_',
                            savepath = './ACPscript/'
                            ):
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form
        
        def ut_to_lt(utctime):
            return self.obsdata.observer.localtime(utctime)
        
        os.makedirs(savepath, exist_ok = True)
        if self.obsdata.schedule == None:
            self.obsdata.schedule = self.obsdata.scheduler(obs_tbl = self.obsdata.target.observable, n_target_for_each_timeslot = n_target_for_each_timeslot, autofocus = autofocus, autofocus_at_init = autofocus_at_init, autofocus_interval_in_hour = autofocus_interval_in_hour, duplicate_when_empty = duplicate_when_empty, delay_minute = delay_minute)

        scheduled_table = self.obsdata.schedule.scheduled
        if len(scheduled_table) > 0:
            # Set observing time 
            time_civil_twilight = self.obsdata.obsnight.sunset_prepare
            time_civil_dawn = self.obsdata.obsnight.sunrise_prepare
            time_astronomical_twilight = self.obsdata.obsnight.sunset_astro
            time_astronomical_dawn = self.obsdata.obsnight.sunrise_astro
            time_obs_start = self.obsdata.obsnight.obs_start
            time_obs_end = self.obsdata.obsnight.obs_end
            
            time_civil_twilight_form = dt_to_form(time_civil_twilight.datetime)
            time_civil_dawn_form = dt_to_form(time_civil_dawn.datetime)
            time_astronomical_twilight_form = dt_to_form(time_astronomical_twilight.datetime)
            time_astronomical_dawn_form = dt_to_form(time_astronomical_dawn.datetime)
            time_obs_start_form = dt_to_form(time_obs_start.datetime)
            time_obs_end_form = dt_to_form(time_obs_end.datetime)
            filename = f'{savepath}{filename_prefix}ACP_{(self.obsdata.obsnight.midnight.datetime).strftime("%y%m%d_%H%M")}_{self.obsdata.name_telescope}.txt'
            
            with open(filename, 'w') as f:
                
                # Preparation
                f.write('; Prepare for the observation (KCT)\n')
                f.write('#WAITUNTIL 1, {}\n'.format(time_civil_twilight_form))
                
                # Cooling
                if ccdcool:
                    f.write('\n; Cooling\n')
                    f.write('#CHILL {}, {}\n'.format(ccdtemp, 1.0))
                
                # Calibration frames
                f.write('\n; Calibration frames\n')
                if bias:
                    f.write("#COUNT {}\n#BIAS\n".format(bias_count)) # bias
                    
                if dark:
                    f.write("#COUNT {}\n#INTERVAL {}\n#DARK\n".format(dark_count, dark_exptime)) # dark
                
                if flat_dusk:
                    f.write('#DUSKFLATS\n') # duskflat
                
                # Image frames
                f.write('\n; Start of the observation\n')
                f.write('#WAITUNTIL 1, {}\n'.format(time_astronomical_twilight_form))
                
                f.write('\n; Targeting\n')
                target_table = Table().from_pandas(scheduled_table.to_pandas())
                target_table = target_table[target_table['id'] != '']
                for target in target_table:
                    if target['obj'] == 'autofocus':
                        f.write('\n#AUTOFOCUS\n')
                    else:
                        name = target['obj']
                        ra = target['ra_hms']
                        dec = target['dec_dms']
                        exptime = target['exptime']
                        filter_ = target['filter']
                        count = target['count']
                        binning = target['binning']
                        f.write('#COUNT {}\n'.format(count))
                        f.write('#INTERVAL {}\n'.format(exptime))
                        f.write('#BINNING {}\n'.format(binning))
                        f.write('#FILTER {}\n'.format(filter_))
                        f.write('#NOPREVIEW \n')
                        f.write('#NOSOLVE\n')
                        f.write('#NOPOINTING\n')
                        f.write('{}\t{}\t{}\n'.format(name,ra,dec))
                        f.write('\n')

                f.write('#QUITAT {}\n'.format(time_astronomical_dawn_form))
                
                if flat_dawn:
                    f.write('#DAWNFLATS\n') # dawnflat
                    
                if shutdown:
                    f.write('\n; Closing\n')
                    f.write('#SHUTDOWN\n')
            print(f'[write_ACPscript_KCT] SAVED! filename : {filename}')
        else:
            print('[write_ACPscript_KCT] FAILED! No observable Target.')

    def show(self, 
             save : bool = False,
             filename_prefix : str = 'ToO_',
             savepath = './ACPscript/'
             ):
        if self.obsdata.schedule == None:
            raise RuntimeError('No specified schedule is defined. Run "write_rts or write_ACPscript"')
        self.obsdata.show(scheduled_table = self.obsdata.schedule.scheduled, save = save, filename_prefix = filename_prefix, savepath = savepath)
    
    def write_log(self,
                  n_target : int = 300,
                  sort_keyword : str  = 'rank',
                  filename_prefix : str = 'ToO_',
                  savepath = './rts/',
                  format_ : str = 'ascii',
                  return_ : bool = True
                  ):
        
        os.makedirs(savepath, exist_ok = True)
        if self.obsdata.schedule == None:
            raise RuntimeError('No specified schedule is defined. Run "write_rts or write_ACPscript" before writing a log file')
        currenttime = self.obsdata.obsnight.midnight.datetime
        all_tbl = self.obsdata.target.all.to_pandas()
        scheduled_tbl = self.obsdata.schedule.all[self.obsdata.schedule.all['scheduled'] == True]
        if len(scheduled_tbl) > 0:
            scheduled_tbl_pd = scheduled_tbl.to_pandas()
            scheduled_tbl = scheduled_tbl_pd.dropna(subset = ['ra'])
            data_table = pd.merge(all_tbl, scheduled_tbl, how = 'left')
            data_table['scheduled'] = np.array(data_table.id.isin(scheduled_tbl.id).tolist())
            log_table = data_table.sort_values(by = ['scheduled',sort_keyword], ascending= [False, True])[:n_target]
            log_table = log_table.sort_values(by = [sort_keyword], ascending = [True])
            log_table = Table().from_pandas(log_table)
        else:
            data_table = all_tbl
            log_table = data_table.sort_values(by = ['scheduled', sort_keyword], ascending= [False, True])[:n_target]
            log_table = log_table.sort_values(by = [sort_keyword], ascending = [True])
            log_table = Table().from_pandas(log_table)
            
        remove_kwargs = ['coord.ra','coord.dec','exec_start','exec_end']
        for remove_kwarg in remove_kwargs:
            if remove_kwarg in log_table.keys():
                log_table.remove_column(remove_kwarg)
            
        filename = f'{savepath}{filename_prefix}{currenttime.strftime("%y%m%d")}_{self.obsdata.name_telescope}.log'
        print(f'[write_log] SAVED! filename : {filename}')
        log_table.write(filename, format = format_, overwrite = True)
        if return_:
            return log_table

################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################

#%% SAMPLE Code (ToO)
tiles = ascii.read('./SkyGridCatalog_KCT_90.csv')
def find_cumulative_index(array, percentile):
    total_sum = sum(array)
    target_cumulative = total_sum * percentile  # 50% of the total sum

    cumulative_sum = 0
    for index, value in enumerate(array):
        cumulative_sum += value
        if cumulative_sum >= target_cumulative:
            return index
    # If the loop completes without finding a suitable index
    return None
idx_50 = find_cumulative_index(tiles['prob_vol_x_stmass'], 0.5)
idx_90 = find_cumulative_index(tiles['prob_vol_x_stmass'], 0.9)
data = tiles[idx_50:idx_90]


if __name__ == '__main__':
    # dirlist = os.listdir('../Archive/')
    # dirlist.sort()
    target = 's230830b'
    name_project = 'GECKO'
    filename_prefix = 'ToO_'
    log_savepath = './log/'
    date = Time.now()
    # ACP config
    ACP_savepath = f'./ACPscript/{target}/'
    # RTS config
    rts_savepath = f'./rts/{target}/'
    n_target_for_each_timeslot = 2

    # def get_isfile_and_data(target, name_telescope):


    #     host_file_key = f'../Archive/{target}/HostGalaxyCatalog_90.csv'
    #     is_host = os.path.isfile(host_file_key)
    #     grid_file_key = f'../Archive/{target}/SkyGridCatalog_{name_telescope}_90.csv'
    #     is_grid = os.path.isfile(grid_file_key)

    #     data_host = Table()
    #     data_targetted = Table()
    #     if is_host:
    #         data_host = ascii.read(host_file_key)
    #     if is_grid:
    #         data_targetted = ascii.read(grid_file_key)
    #     result = dict()
    #     result['host'] = dict()
    #     result['host']['exist'] = is_host
    #     result  = data_host
    #     result['grid'] = dict()
    #     result['grid']['exist'] = is_grid
    #     result  = data_targetted
    #     return result

    # print(f'len(targetted) = {len(get_isfile_and_data(target=target, name_telescope="KCT")["host"]["data"])}')
    
    
    ####################### 
    # KCT (ACP)
    #######################
    
    name_telescope = 'KCT'
    data = ascii.read('./SkyGridCatalog_KCT_90_select.csv')
    data['filter'] = "r"
    data['count'] = '6'
    data['exptime'] = '300'
    date = Time.now()
    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """

    scheduler_host = ObsScheduler(target_db= data,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_ACPscript_KCT(filename_prefix= filename_prefix, savepath = ACP_savepath, shutdown = False)
    scriptmaker_host.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)

    scheduler_grid = ObsScheduler(target_db= data,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_grid = ScriptMaker(scheduler_grid)
    # Action
    scriptmaker_grid.write_ACPscript_KCT(filename_prefix= filename_prefix, savepath = ACP_savepath, shutdown = False)
    scriptmaker_grid.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_grid.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)
    # See log file to check the observability of the targets
    plt.figure(dpi = 300)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
    plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    #plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
    #plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    maintarget = mainTarget(name_telescope=scheduler_host.name_telescope, name_project= scheduler_host.name_project, observer = scheduler_host.observer, target_ra = scheduler_host.target.all['ra'][0], target_dec = scheduler_host.target.all['dec'][0])
    maintarget.staralt(utctime = date)
    
    #######################
    # RASA36 (ACP)
    #######################
    
    name_telescope = 'RASA36'
    data = ascii.read('./test/SkyGridCatalog_RASA36_90.csv')
    #data = ascii.read('./SkyGridCatalog_RASA36_90_fromKMTNet.csv')
    #data['count'] = '30'
    #data['exptime'] = '60'
    date = Time.now()

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    scheduler_host = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_ACPscript_RASA36(filename_prefix= filename_prefix, savepath = ACP_savepath)
    rasalog = scriptmaker_host.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = True)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)

    scheduler_grid = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_grid = ScriptMaker(scheduler_grid)
    # Action
    scriptmaker_grid.write_ACPscript_RASA36(filename_prefix= filename_prefix, savepath = ACP_savepath)
    scriptmaker_grid.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_grid.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)
    # See log file to check the observability of the targets
    # Check
    plt.figure(dpi = 300)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
    plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    #plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
    #plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()

    
    #######################
    # LSGT (ACP)
    #######################
    name_telescope = 'LSGT'
    data = ascii.read('./HostGalaxyCatalog_90.csv')

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    scheduler_host = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_ACPscript_LSGT(filename_prefix= filename_prefix, savepath = ACP_savepath, period_script= 3)
    scriptmaker_host.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)

    scheduler_grid = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_grid = ScriptMaker(scheduler_grid)
    # Action
    scriptmaker_grid.write_ACPscript_LSGT(filename_prefix= filename_prefix, savepath = ACP_savepath, period_script= 3)
    scriptmaker_grid.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_grid.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)
    # See log file to check the observability of the targets
    # Check
    plt.figure(dpi = 300)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
    plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
    plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    maintarget = mainTarget(name_telescope=scheduler_host.name_telescope, name_project= scheduler_host.name_project, observer = scheduler_host.observer, target_ra = scheduler_host.target.all['ra'][0], target_dec = scheduler_host.target.all['dec'][0])
    maintarget.staralt(utctime = date)
    
    
    ####################### 
    # LOAO (2x RTS)
    #######################
    
    name_telescope = 'LOAO'
    data = ascii.read('./HostGalaxyCatalog_90.csv')

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    scheduler_host = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_host.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)

    scheduler_grid = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_grid = ScriptMaker(scheduler_grid)
    # Action
    scriptmaker_grid.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_grid.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_grid.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)
    # See log file to check the observability of the targets
    # Check
    plt.figure(dpi = 300)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
    plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
    plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    maintarget = mainTarget(name_telescope=scheduler_host.name_telescope, name_project= scheduler_host.name_project, observer = scheduler_host.observer, target_ra = scheduler_host.target.all['ra'][0], target_dec = scheduler_host.target.all['dec'][0])
    maintarget.staralt(utctime = date)

    ####################### 
    # CBNUO (2x RTS)
    #######################
    
    name_telescope = 'CBNUO'
    data = ascii.read('./S230521k/SkyGridCatalog_CBNUO_90.csv')

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    scheduler_host = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_host.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)

    scheduler_grid = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_grid = ScriptMaker(scheduler_grid)
    # Action
    scriptmaker_grid.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_grid.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_grid.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)
    # See log file to check the observability of the targets
    # Check
    plt.figure(dpi = 300)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
    plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
    plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    maintarget = mainTarget(name_telescope=scheduler_host.name_telescope, name_project= scheduler_host.name_project, observer = scheduler_host.observer, target_ra = scheduler_host.target.all['ra'][0], target_dec = scheduler_host.target.all['dec'][0])
    maintarget.staralt(utctime = date)
    
    
    ####################### 
    # SAO (2x RTS)
    #######################
    
    name_telescope = 'SAO'
    data = ascii.read('./SkyGridCatalog_KMTNet_90.csv')

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    scheduler_host = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_host.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)

    scheduler_grid = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_grid = ScriptMaker(scheduler_grid)
    # Action
    scriptmaker_grid.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_grid.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_grid.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)
    # See log file to check the observability of the targets
    # Check
    plt.figure(dpi = 300)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
    plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
    plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    maintarget = mainTarget(name_telescope=scheduler_host.name_telescope, name_project= scheduler_host.name_project, observer = scheduler_host.observer, target_ra = scheduler_host.target.all['ra'][0], target_dec = scheduler_host.target.all['dec'][0])
    maintarget.staralt(utctime = date)
    
    
    #######################
    # CTIO (1x RTS)
    #######################
    
    name_telescope = 'KMTNet'
    data = ascii.read('./SkyGridCatalog_KMTNet_90.csv')
    data = ascii.read('./rts/MS181101ab_INITIAL/230519_KMTNet_CTIO.log')

    
    n_target_for_each_timeslot = 1 ##### IMPORTANT?? #####
    name_telescope = 'KMTNet_CTIO'
    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """

    scheduler_host = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_host.write_log(n_target = 10000, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)


    scheduler_grid = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_grid = ScriptMaker(scheduler_grid)
    # Action
    scriptmaker_grid.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_grid.write_log(n_target = 10000, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_grid.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)
    # See log file to check the observability of the targets
    # Check
    plt.figure(dpi = 300)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
    plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
    plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    maintarget = mainTarget(name_telescope=scheduler_host.name_telescope, name_project= scheduler_host.name_project, observer = scheduler_host.observer, target_ra = scheduler_host.target.all['ra'][0], target_dec = scheduler_host.target.all['dec'][0])
    maintarget.staralt(utctime = date)


    ####################### 
    # SAAO (1x RTS)
    #######################
    
    name_telescope = 'KMTNet'
    data = ascii.read('./SkyGridCatalog_KMTNet_90.csv')
    n_target_for_each_timeslot = 1 ##### IMPORTANT?? #####
    name_telescope = 'KMTNet_SAAO'
    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """

    scheduler_host = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_host.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)


    scheduler_grid = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_grid = ScriptMaker(scheduler_grid)
    # Action
    scriptmaker_grid.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_grid.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_grid.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)
    # See log file to check the observability of the targets
    # Check
    plt.figure(dpi = 300)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
    plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
    plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    maintarget = mainTarget(name_telescope=scheduler_host.name_telescope, name_project= scheduler_host.name_project, observer = scheduler_host.observer, target_ra = scheduler_host.target.all['ra'][0], target_dec = scheduler_host.target.all['dec'][0])
    maintarget.staralt(utctime = date)
    
    
    ######################## 
    # SSO (1x RTS)
    #######################
    
    name_telescope = 'KMTNet'
    data = ascii.read('./SkyGridCatalog_KMTNet_90.csv')

    n_target_for_each_timeslot = 1 ##### IMPORTANT?? #####
    name_telescope = 'KMTNet_SSO'
    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    scheduler_host = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_host.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)

    scheduler_grid = ObsScheduler(target_db= data ,
                                    date = date,
                                    name_project = name_project,
                                    name_telescope = name_telescope,
                                    entire_night = False)
    # Define target
    scriptmaker_grid = ScriptMaker(scheduler_grid)
    # Action
    scriptmaker_grid.write_rts(filename_prefix= filename_prefix, savepath = rts_savepath, n_target_for_each_timeslot= n_target_for_each_timeslot)
    scriptmaker_grid.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= rts_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_grid.show(save = True, filename_prefix = filename_prefix, savepath = rts_savepath)
    # See log file to check the observability of the targets
    # Check
    plt.figure(dpi = 300)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
    plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
    plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    maintarget = mainTarget(name_telescope=scheduler_host.name_telescope, name_project= scheduler_host.name_project, observer = scheduler_host.observer, target_ra = scheduler_host.target.all['ra'][0], target_dec = scheduler_host.target.all['dec'][0])
    maintarget.staralt(utctime = date)


#%% SAMPLE (IMSNG)
if __name__ == '__main__':
    ####################### KCT (IMSNG)
    n_date = 100
    for i in range(n_date):
        date = Time.now() + i *u.day
        ACP_savepath = f'./IMSNG/'
        name_telescope = 'KCT'
        name_project = 'IMSNG'
        filename_prefix = 'IMSNG_'
        duplicate_when_empty = True
        data = ascii.read('./alltarget_prior15.dat', format = 'fixed_width')
        data['weight'] = data['priority']

        scheduler_host = ObsScheduler(target_db= data,
                                            date = date,
                                            name_project = name_project,
                                            name_telescope = name_telescope,
                                            entire_night = True)
        # Define target
        scriptmaker_host = ScriptMaker(scheduler_host)
        # Action
        scriptmaker_host.write_ACPscript_KCT(filename_prefix= filename_prefix, savepath = ACP_savepath, shutdown = True, duplicate_when_empty= duplicate_when_empty, bias= True, dark = True)
        #scriptmaker_host.write_log(n_target = 300, sort_keyword = 'priority', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = False)
        #scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)
#%%

    ######################## LSGT (IMSNG)
    date = Time.now() 
    ACP_savepath = f'./IMSNG/'
    name_telescope = 'LSGT'
    project = 'IMSNG'
    filename_prefix = 'IMSNG_'
    duplicate_when_empty = True
    data = ascii.read('./alltarget_prior2.dat', format = 'fixed_width')

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    scheduler_host = ObsScheduler(target_db= data,
                                        date = date,
                                        name_project = name_project,
                                        name_telescope = name_telescope,
                                        entire_night = True)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_ACPscript_LSGT(filename_prefix= filename_prefix, savepath = ACP_savepath, duplicate_when_empty= duplicate_when_empty)
    scriptmaker_host.write_log(n_target = 300, sort_keyword = 'priority', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)
    
    ######################## RASA36 (IMSNG)
    date = Time.now()
    ACP_savepath = f'./IMSNG/'
    name_telescope = 'RASA36'
    project = 'IMSNG'
    filename_prefix = 'IMSNG_'
    duplicate_when_empty = True
    data = ascii.read('./alltarget_prior2.dat', format = 'fixed_width')

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    scheduler_host = ObsScheduler(target_db= data,
                                        date = date,
                                        name_project = name_project,
                                        name_telescope = name_telescope,
                                        entire_night = True)
    # Define target
    scriptmaker_host = ScriptMaker(scheduler_host)
    # Action
    scriptmaker_host.write_ACPscript_RASA36(filename_prefix= filename_prefix, savepath = ACP_savepath, shutdown = True, duplicate_when_empty= duplicate_when_empty)
    scriptmaker_host.write_log(n_target = 300, sort_keyword = 'priority', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = False)
    scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)
