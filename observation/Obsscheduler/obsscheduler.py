


#%%
import os
import sys

# 현재 스크립트 파일의 디렉토리 경로
dir_path = os.path.dirname(os.path.realpath(__file__))

# 'a'와 'b' 모듈을 import하기 위해 상위 폴더 경로를 추가
sys.path.append(os.path.join(dir_path, '..'))

from mainobserver import mainObserver
from mainconfig import mainConfig
from maintarget import mainTarget
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astroplan import *
# from astroplan import AltitudeConstraint, AirmassConstraint, MoonSeparationConstraint, AtNightConstraint
# from astroplan import observability_table, FixedTarget
from astropy.time import Time
import datetime
import uuid
import matplotlib.pyplot as plt
from astropy.table import join
import os
# %%
class ObsScheduler(mainConfig):
    """
    An observation scheduler class that generates observing schedules based on a target database, time constraints, and the telescope name.
    - History
    (23.04.22) Written by Hyeonho Choi 
    (23.04.23) Change in write_ACPscript_RASA36, write_ACPscript_KCT : change the order to apply #NOSOLVING for all targets
    
    Parameters
    ==========
    1. target_db : Astropy.table.Table, pandas.Dataframe, Dict
        data containing target information.
        - Must indluded keys
        ========
        obj : object name
        ra : target right ascension, '10:00:00', "10 00 00", "10h 00m 00s", '150.0000'
        dec : target declination, '-30:00:00', "-30 00 00", "-30d 00m 00s", '-30.0000'
        Other keys (filters, exptime, ...) will be set default(from Scheduler_telescope.config) if not provided
        ========
        
    2. date : datetime.datetime, astropy.time.Time, str (default=Time.now())
        The observation date (in UT).
        - Current supported string format
        ========
        str : 20230421, 2023/04/21, 2023.04.21, 230421, 23/04/21, 23.04.21
        astropy.time.Time
        datetime.datetime
        ========
        
    3. name_telescope : str (default='KCT')
        The name of the telescope for which the observation plan is being made.
        - Current supported telescopes
        ========
        name_telescope : 'KCT', 'RASA36', 'CBNUO', 'LSGT'
        ========
        
    4. entire_night : bool (default=False)
        Whether to generate a plan for the entire night or only part of it.
        
    5. duplicate_when_empty : bool (default=True)
        Whether to duplicate the observation plan when there are no scheduled targets.

    Methods
    =======
    1. scorer(obs_tbl: Table, telescope_altaz: SkyCoord, utctime: Time, datetime(Time.now()), duplicate: bool(False)) -> np.array
        Calculates scores for the given observation table based on various constraints and weights.
    2. scheduler(delay_time : float) -> Table
        Creates a schedule for the observation based on the observable targets and constraints.
    3. show(save: bool, filename_prefix: str, savepath: bool) -> None
        Visualize observation schedule
    4. write_Gsheet(spreadsheet_url: str, authorize_jsonfile: str, sheetname_prefix: str) -> None
        Writes the observation plan to google spreadsheet
    5. write_rts(scheduled_only: bool, filename_prefix: str, savepath: str) -> None
        Writes the observation rts
    6. write_ACPscript_KCT(filename_prefix: str, savepath: str, shutdown : bool) -> None
        Writes the observation plan in ACP script format for KCT telescope.
    7. write_ACPscript_RASA36(filename_prefix: str, savepath: str) -> None
        Writes the observation plan in ACP script format for RASA36 telescope.
    8. write_ACPscript_LSGT(filename_prefix: str, savepath: str, period_script: int) -> None
        Writes the observation plan in ACP script format for LSGT telescope.
    9. write_NINAscript() : TBD
    """

    def __init__(self,
                 target_db,
                 date = Time.now(),
                 name_telescope = 'KCT',
                 entire_night : bool = False,
                 duplicate_when_empty : bool = True
                 ):

        super().__init__(name_telescope = name_telescope)
        self.constraints = self._set_constraints(**self.config)
        self.score_weight = self._set_weight_score(**self.config)
        
        self._name_telescope = name_telescope
        self._target_db = target_db
        self.all_target_tbl = self._match_targetDB_format()
        self._date = date
        self._entire_night = entire_night
        self._duplicate = duplicate_when_empty
        
        self.observer_tcspy = mainObserver(name_telescope= name_telescope)
        self.observer_astroplan = self.observer_tcspy._observer
        self.date = self._match_date_format()
        self.tonight = self._get_tonight(horizon = -18)
        self.midnight = self._get_midnight()
        self.civil_night = self._get_tonight(horizon = 0)
        self.moon_radec = self._get_moon_radec(self.date)
        self.obs_tbl = self._get_observable_target()
        self.obs_tbl = self._update_tbl(self.obs_tbl)
        self.scheduled_table = None
    
    def scorer(self,
               obs_tbl : Table,
               telescope_altaz : SkyCoord,
               utctime : datetime or Time = None,
               duplicate : bool = False):
        """
        Calculates scores for the given observation table based on various constraints and weights.

        Parameters
        ==========
        1. obs_tbl : Table
            The observation table containing target-related information.
            required columns : ra, dec, priority, transit, maxalt, scheduled, 
        2. telescope_altaz : SkyCoord
            The current AltAz coordinates of the telescope.
        3. utctime : datetime or Time, optional, default=None
            The time to be used for the calculations. If None, the current time will be used.
        4. duplicate : bool, optional, default=False
            If True, duplicate observations will be considered in the scoring.

        Returns
        =======
        1. score : ndarray
            An array of scores for each target in the observation table.
        """
        
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)

        all_target_radec = self._match_coord_format(ra = obs_tbl['ra'], dec = obs_tbl['dec'])
        all_target_altaz = self.observer_astroplan.altaz(utctime, target = all_target_radec)
        all_target_alt = all_target_altaz.alt.value
        all_target_priority = obs_tbl['priority'].astype(float)
        all_target_hourangle_tmp = utctime.datetime - obs_tbl['transit']
        all_target_hourangle = np.array([target_hourangle.total_seconds() for target_hourangle in all_target_hourangle_tmp])
        target_sepdist = SkyCoord.separation(telescope_altaz, all_target_altaz).value
        
        moon_altaz = self.observer_astroplan.altaz(utctime, target = self.moon_radec)
        moonsep = SkyCoord.separation(all_target_altaz, moon_altaz).value
        
        score = np.ones(len(obs_tbl))
        # constraints
        constraint_moonsep = moonsep > self.config['TARGET_MOONSEP']
        constraint_altitude = all_target_alt > self.config['TARGET_MINALT'] +5
        constraint_meridian = (np.abs(all_target_hourangle) > self.config['TARGET_MERIDIAN_SEC']) #(-self.config['TARGET_MERIDIAN_SEC'] > all_target_hourangle) |  ( all_target_hourangle > 0)
        constraint_maxslew = target_sepdist < self.config['TARGET_MAXSLEW']
        constrint_duplicate = ~obs_tbl['scheduled']
        score *= constraint_moonsep
        score *= constraint_altitude
        score *= constraint_meridian
        score *= constraint_maxslew
        if not duplicate:
            score *= constrint_duplicate    
        # score
        score_relative_alt = self.score_weight['TARGET_WEIGHT_RELATIVE_ALT'] * all_target_alt / obs_tbl['maxalt']
        score_absolute_alt = self.score_weight['TARGET_WEIGHT_ABSOLUTE_ALT'] * all_target_alt / 90
        score_priority = self.score_weight['TARGET_WEIGHT_PRIORITY']* (1- all_target_priority / np.max(all_target_priority))
        score_slewtime = self.score_weight['TARGET_WEIGHT_SLEW'] * (target_sepdist / self.config['OVERHEAD_SPEED_SLEWING']) / np.max((target_sepdist / self.config['OVERHEAD_SPEED_SLEWING']))
        score_all = (score_relative_alt + score_absolute_alt + score_priority + score_slewtime) / 4# (self.score_weight['TARGET_WEIGHT_RELATIVE_ALT'] + self.score_weight['TARGET_WEIGHT_ABSOLUTE_ALT'] + self.score_weight['TARGET_WEIGHT_PRIORITY'] + self.score_weight['TARGET_WEIGHT_SLEW'])
        score *= score_all
        return score
    
    def scheduler(self,
                  delay_time : float = 0
                  ):
        """
        Schedules the observations based on provided constraints, delays, and priorities.

        Parameters
        ==========
        1. delay_time : float, optional, default=0
            The time delay in minutes for the start of the observation.

        Returns
        =======
        1. scheduled_table : Table
            A table containing the scheduled observations and their details.
        """
        
        # Set observing time 
        time_astronomical_dawn = self.observer_tcspy.sun_risetime(self.date, mode = 'next', horizon= -18)
        time_astronomical_twilight = self.observer_tcspy.sun_settime(utctimes = time_astronomical_dawn, mode = 'previous', horizon= -18)
        time_civil_dawn = self.observer_tcspy.sun_risetime(time_astronomical_dawn, mode = 'next', horizon= 0)
        time_civil_twilight = self.observer_tcspy.sun_settime(utctimes = time_civil_dawn, mode = 'previous', horizon= 0)
        
        #time_astronomical_twilight = self.tonight[0] + delay_time * u.minute
        time_obs_start = self.tonight[0] + delay_time * u.minute
        time_obs_end = self.tonight[1] - delay_time * u.minute
        # Set transitioner
        target_bias = self._transitioner(name = 'bias', exptime = self.config['OVERHEAD_TIME_READOUT'], count = self.config['SCHEDULER_BIASCOUNT'], id_ ='')
        target_dark = self._transitioner(name = 'dark', exptime = self.config['SCHEDULER_DARKEXPTIME'], count = self.config['SCHEDULER_DARKCOUNT'], id_ ='')
        target_flat_dusk = self._transitioner(name = 'duskflat', exptime = self.config['SCHEDULER_FLATEXPTIME'], count = self.config['SCHEDULER_FLATCOUNT'], filter = self.config['SCHEDULER_FLATFILTSET'], id_ ='')
        target_flat_dawn = self._transitioner(name = 'dawnflat', exptime = self.config['SCHEDULER_FLATEXPTIME'], count = self.config['SCHEDULER_FLATCOUNT'], filter = self.config['SCHEDULER_FLATFILTSET'], id_ ='')
        
        obs_tbl = self.obs_tbl.copy()
        scheduled_table = Table(names = ['obj', 'ra', 'dec', 'ra_deg', 'dec_deg', 'time_start', 'time_end', 'exptime', 'filter', 'count', 'binning', 'priority', 'score', 'id'], 
                                dtype = ['object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object'])
                                #dtype = ['str', 'object', 'object', 'float', 'float', 'object', 'object', 'str', 'str', 'str', 'str', 'float', 'float', 'str'])
        
        # Preparation of the observation
        time_observation = time_civil_twilight
        if time_obs_start < time_civil_twilight + 30 * u.minute:
            # Dusk flat
            if self.config['SCHEDULER_FLAT_DUSK']:
                scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_flat_dusk, score = 1)
                target_overhead_flat = self._transitioner(name = f'overhead_flat', exptime = self._calculate_tot_overhead(target_flat_dusk), count = 1)
                scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_overhead_flat)
            # Bias
            if self.config['SCHEDULER_BIAS']:
                scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_bias, score = 1)
                target_overhead_bias = self._transitioner(name = f'overhead_bias', exptime = self._calculate_tot_overhead(target_bias), count = 1)
                scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_overhead_bias)
            # Dark
            if self.config['SCHEDULER_DARK']:
                scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_dark, score = 1)
                target_overhead_dark = self._transitioner(name = f'overhead_dark', exptime = self._calculate_tot_overhead(target_dark), count = 1)
                scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_overhead_dark)
            
        time_observation = time_obs_start
        autofocus_time = time_observation - 5 * u.hour
        telescope_altaz = obs_tbl['target_tcspy'][0].altaz(time_observation)
        target_best = None
        while time_observation < time_obs_end:
            # Autofocus
            if self.config['SCHEDULER_AUTOFOCUS']:
                if (time_observation - autofocus_time > self.config['SCHEDULER_AFTIME'] * u.hour) * (time_astronomical_dawn - time_observation > self.config['SCHEDULER_AFTIME'] * u.hour):
                    target_autofocus = self._transitioner(name = 'autofocus', exptime = self.config['OVERHEAD_TIME_AUTOFOCUS'], count = 1)
                    scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_autofocus, score = 1)
                    autofocus_time = time_observation
            # Target selection 
            score = self.scorer(obs_tbl, telescope_altaz, utctime = time_observation)

            if np.sum(score) == 0:
                
                if self._duplicate:
                    if not target_best is None:
                        obs_tbl_all = self.obs_tbl.copy()
                        obs_tbl_all = obs_tbl_all[~(target_best['obj'] == obs_tbl_all['obj'])] # Prevent subsequent observations 
                    score = self.scorer(obs_tbl_all, telescope_altaz, utctime = time_observation, duplicate= True)
                    idx_best = np.argmax(score)
                    target_best = obs_tbl_all[idx_best]
                    target_score = score[idx_best]
                else:
                    target_best = self._transitioner(name = f'empty_target', exptime = 600, count = 1)
                    target_score = -999
            else:
                idx_best = np.argmax(score)
                target_best = obs_tbl[idx_best]
                target_score = score[idx_best]
                
            scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_best, score = target_score)
            target_overhead = self._transitioner(name = f'overhead_{target_best["obj"]}', exptime = self._calculate_tot_overhead(target_best), count = 1)
            scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_overhead)
            obs_tbl = self._update_status_obs_tbl(obs_tbl = obs_tbl, target = target_best)
        
        if self.config['SCHEDULER_FLAT_DAWN']:
            scheduled_table, time_observation = self._insert_schedule(start_time = time_astronomical_dawn, schedule_table=scheduled_table, target = target_flat_dawn, score = 1)
            target_overhead_flat = self._transitioner(name = f'overhead_flat', exptime = self._calculate_tot_overhead(target_flat_dawn), count = 1)
            scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_overhead_flat)
        scheduled_table = scheduled_table[scheduled_table['score'] > 0 ]
        self.scheduled_table = scheduled_table
        return scheduled_table
    
    def show(self,
             save : bool = False,
             filename_prefix : str = '',
             savepath = './ACPscript/'
             ):
        """
        Displays and saves the observing plan as a graph with altitude over time.

        Parameters
        ==========
        1. save : bool, optional, default=False
            Whether to save the plot to a file.
        2. filename_prefix : str, optional, default=''
            A prefix for the filename of the saved plot.
        3. savepath : str, optional, default='./'
            The directory where the plot will be saved.

        Returns
        =======
        None, but it displays the plot and saves it to a file if 'save' is set to True.
        """

        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        n_obj = len(set(scheduled_table['obj']))
        target_table = scheduled_table[scheduled_table['ra'] != '0000']
        time_astronomical_dawn = self.observer_tcspy.sun_risetime(self.date, mode = 'next', horizon= -18)
        time_astronomical_twilight = self.observer_tcspy.sun_settime(utctimes = time_astronomical_dawn, mode = 'previous', horizon= -18)
        time_civil_dawn = self.observer_tcspy.sun_risetime(time_astronomical_dawn, mode = 'next', horizon= 0)
        time_civil_twilight = self.observer_tcspy.sun_settime(utctimes = time_civil_dawn, mode = 'previous', horizon= 0)

        time_range_start, time_range_end = time_civil_twilight.datetime - datetime.timedelta(hours = 2), time_civil_dawn.datetime + datetime.timedelta(hours = 2)
        time_axis = np.arange(time_range_start, time_range_end, datetime.timedelta(minutes = 5))
        moon_altaz = self.observer_tcspy.moon_altaz(time_axis)
        sun_altaz = self.observer_tcspy.sun_altaz(time_axis)
        
        plt.figure(dpi = 300, figsize = (10, 4))
        plt.title(f'Observing plan[{self.name_telescope}]')
        for target in target_table:
            target_id = target['id']
            target_tcspy = self.obs_tbl[self.obs_tbl['id'] == target_id]['target_tcspy'][0]
            target_altaz = target_tcspy.altaz(time_axis)
            obs_time_axis = np.arange(target['time_start'], target['time_end'], datetime.timedelta(minutes = 1))
            obs_target_altaz = target_tcspy.altaz(obs_time_axis)
            plt.plot(target_altaz.obstime.datetime, target_altaz.alt.value, c='k', linewidth = 0.3)
            plt.plot(obs_target_altaz.obstime.datetime, obs_target_altaz.alt.value, c='r')
        if (self.date.datetime < time_range_end + datetime.timedelta(hours = 3)) & (self.date.datetime > time_range_start - datetime.timedelta(hours = 3)) * (~self._entire_night):
            plt.axvline(self.date.datetime, linestyle = '--', c='r', label = 'Now')
        plt.plot(moon_altaz.obstime.datetime, moon_altaz.alt.value, c = 'b', label ='Moon', linestyle = ':')
        plt.plot(sun_altaz.obstime.datetime, sun_altaz.alt.value, c = 'r', label = 'Sun', linestyle = ':')
        
        plt.fill_betweenx([10,90], time_astronomical_twilight.datetime, time_astronomical_dawn.datetime, alpha = 0.1, color='k')
        plt.fill_betweenx([10,90], time_civil_twilight.datetime, time_civil_dawn.datetime, alpha = 0.05, color='k')
        plt.axvline(x=time_astronomical_twilight.datetime, linestyle = '-', c='k', linewidth = 0.5)
        plt.axvline(x=time_astronomical_dawn.datetime, linestyle = '-', c='k', linewidth = 0.5)
        plt.axvline(x=time_civil_dawn.datetime, linestyle = '--', c='k', linewidth = 0.5)
        plt.axvline(x=time_civil_twilight.datetime, linestyle = '--', c='k', linewidth = 0.5)
        plt.axhline(y=self.config['TARGET_MINALT'], linestyle = '--', c='r', linewidth = 1)
        plt.text(time_civil_twilight.datetime-datetime.timedelta(minutes=140), 85, f'N_target = {len(scheduled_table)}[{n_obj}]', fontsize = 10)
        plt.text(time_astronomical_twilight.datetime, 92, 'Ast.S.set', fontsize = 10)
        plt.text(time_civil_twilight.datetime-datetime.timedelta(minutes=00), 92, 'S.set', fontsize = 10)
        plt.text(time_civil_dawn.datetime-datetime.timedelta(minutes=00), 92, 'S.rise', fontsize = 10)
        plt.xlim(time_range_start - datetime.timedelta(hours = 1), time_range_end + datetime.timedelta(hours = 1))
        plt.ylim(10, 90)
        plt.legend(loc = 1)
        plt.xlabel('UTC [mm-dd hh]')
        plt.ylabel('Altitude [deg]')
        plt.grid()
        plt.xticks(rotation = 45)
        if save:
            os.makedirs(savepath, exist_ok = True)
            filename = f'{savepath}{filename_prefix}ObservingSchedule_{(self.midnight.datetime).strftime("%y%m%d")}_{self.name_telescope}.png'
            plt.savefig(filename, bbox_inches = 'tight')
            
    def write_Gsheet(self,
                     spreadsheet_url : str,
                     authorize_jsonfile : str,
                     sheetname_prefix : str = '',
                     ):
        """
        Writes the observation plan to a Google Sheet.

        Parameters
        ==========
        1. spreadsheet_url : str
            The URL of the Google Sheets document.
        2. authorize_jsonfile : str
            The path to the JSON file containing the authorization credentials.
        3. sheetname_prefix : str, optional, default=''
            A prefix to add to the sheet name when creating a new sheet in the Google Sheets document.

        Returns
        =======
        None
        """

        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        # To use this, a token is required for authorization for your google account! 
        from googlesheet import GoogleSheet
        scheduled_table['time_start'] = [str(target['time_start']) for target in scheduled_table]
        scheduled_table['time_end'] = [str(target['time_end']) for target in scheduled_table]
        gs = GoogleSheet(spreadsheet_url= spreadsheet_url, authorize_json_file= authorize_jsonfile )
        print(gs.address_sheets)
        sheetname = f'{sheetname_prefix}{self.midnight.datetime.strftime("%y%m%d")}_{self.name_telescope}'
        gs.write_sheet_data(sheet_name = sheetname, data = scheduled_table, append = False )
    
    def write_rts(self,
                  scheduled_only : bool = True,
                  filename_prefix : str = '',
                  savepath = './rts/'
                  ):
        """
        Writes the observation plan to a text file.

        Parameters
        ==========
        1. scheduled_only : bool, optional, default=True
            If True, only scheduled targets will be written. If False, all targets will be written.
        2. filename_prefix : str, optional, default=''
            A prefix to add to the filename when creating a new file.
        3. savepath : str, optional, default='./'
            The path to the directory where the file will be saved.

        Returns
        =======
        None
        """
        
        os.makedirs(savepath, exist_ok = True)
        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        def ut_to_lt(utctime):
            return self.observer_tcspy.localtime(utctime)
        def dt_to_form_hoursec(dt):
            form ='{:02}:{:02}'.format(dt.hour, dt.minute)
            return form
        if scheduled_only:
            obs_info = join(scheduled_table, self.obs_tbl)
        else:
            obs_info = self.obs_tbl
        obs_info.sort('ra_deg')
        # Target
        all_target_radec = self._match_coord_format(ra = obs_info['ra'], dec = obs_info['dec'])
        all_target_altaz = self.observer_astroplan.altaz(self.midnight, target = all_target_radec)
        risetime_lt = [dt_to_form_hoursec(ut_to_lt(utctime)) for utctime in obs_info['rise']]
        transittime_lt = [dt_to_form_hoursec(ut_to_lt(utctime)) for utctime in obs_info['transit']]
        settime_lt = [dt_to_form_hoursec(ut_to_lt(utctime)) for utctime in obs_info['set']]
        moon_altaz = self.observer_astroplan.altaz(self.midnight, target = self.moon_radec)
        moonsep = SkyCoord.separation(all_target_altaz, moon_altaz).value.astype(int)
        # Moon
        moon_radec = self.observer_tcspy.moon_radec(self.midnight)
        moonphase = round(self.observer_tcspy.moon_phase(self.midnight),2)
        moon_ra_str = f'{int(moon_radec.ra.hms.h)}:{int(moon_radec.ra.hms.m)}:{round(moon_radec.ra.hms.s,2)}'
        moon_dec_str = f'{int(moon_radec.dec.dms.d)}:{np.abs(int(moon_radec.dec.dms.m))}:{round(np.abs(moon_radec.dec.dms.s),2)}'
        moon_set_sepration = int(self.config['TARGET_MOONSEP'])
        # Sun   
        sunrisetime = ut_to_lt(self.observer_tcspy.sun_risetime(self.date, mode = 'next', horizon= -18).datetime)
        sunsettime = ut_to_lt(self.observer_tcspy.sun_settime(utctimes = sunrisetime, mode = 'previous', horizon= -18).datetime)
        tonight = ut_to_lt(self.midnight.datetime).strftime('%Y-%m-%d')
        filename = f'{savepath}{filename_prefix}rts_{ut_to_lt(self.midnight.datetime).strftime("%y%m%d")}_{self.name_telescope}.txt'

        # Target table
        rts_table = Table()
        rts_table['name'] = obs_info['obj']
        rts_table['ra'] = obs_info['ra']
        rts_table['dec'] = obs_info['dec']
        rts_table['rise(LT)'] = risetime_lt
        rts_table['transit(LT)'] = transittime_lt
        rts_table['set(LT)'] = settime_lt
        rts_table['moon_dist(deg)'] = moonsep
        rts_table['priority'] = obs_info['priority']
        #rts_table['filter'] = obs_info['filter']
        #rts_table['exptime'] = obs_info['exptime']
        #rts_table['count'] = obs_info['count']
        #rts_table['binning'] = obs_info['binning']
        # Sorting
        #rts_table.sort('priority')
        rts_table.write(filename, format = 'ascii.fixed_width', overwrite = True)
        
        # add info
        f = open(filename, 'r')
        original_data = f.readlines()
        with open(filename, "w") as f:
            f.write('#\tObserving Date = '+tonight+'\n')
            f.write('#\tObservatory = '+self.name_telescope+'\n')
            f.write('#\t-18 deg sunset = '+dt_to_form_hoursec(sunsettime)+'\n')
            f.write('#\t-18 deg sunrise = '+dt_to_form_hoursec(sunrisetime)+'\n')
            f.write('#\tMoon ra, dec = '+ moon_ra_str +'   '+ moon_dec_str+'\n')
            f.write('#\tMoon phase = '+str(moonphase)+'\n')
            f.write('#\tMoon seperation limit = '+str(moon_set_sepration)+'\n')
            f.write('#\tAltitude limit = '+str(self.config['TARGET_MINALT'])+'\n')
            f.write('='*108+'\n')
            for line in original_data:
                f.write(line)
            f.write('='*108+'\n')

    def write_ACPscript_KCT(self,
                            filename_prefix : str = '',
                            savepath = './ACPscript/',
                            shutdown : bool = True,
                            ):
        """
        Writes the observation plan in ACP script format for KCT (Korea Cerro Tololo) telescope.

        Parameters
        ==========
        1. filename_prefix : str, optional, default=''
            A prefix to add to the filename when creating a new file.
        2. savepath : str, optional, default='./'
            The path to the directory where the file will be saved.
        3. shutdown : bool, optional, default=True
            If True, the script will include a shutdown command at the end.

        Returns
        =======
        None
        """
        
        os.makedirs(savepath, exist_ok = True)
        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form
        
        # Set observing time 
        time_civil_twilight = self.civil_night[0] 
        time_civil_dawn = self.civil_night[1] 
        time_astronomical_twilight = self.tonight[0]
        time_astronomical_dawn = self.tonight[1]
        
        time_civil_twilight_form = dt_to_form(time_civil_twilight.datetime)
        time_civil_dawn_form = dt_to_form(time_civil_dawn.datetime)
        time_astronomical_twilight_form = dt_to_form(time_astronomical_twilight.datetime)
        time_astronomical_dawn_form = dt_to_form(time_astronomical_dawn.datetime)
        filename = f'{savepath}{filename_prefix}ACP_{(self.midnight.datetime).strftime("%y%m%d_%H%m")}_{self.name_telescope}.txt'
        
        with open(filename, 'w') as f:
            
            # Preparation
            f.write('; Prepare for the observation\n')
            f.write('#WAITUNTIL 1, {}\n'.format(time_civil_twilight_form))
            
            # Cooling
            if self.config['SCHEDULER_CCDCOOL']:
                f.write('\n; Cooling\n')
                f.write('#CHILL {}, {}\n'.format(self.config['SCHEDULER_CCDTEMP'], 1.0))
            
            # Calibration frames
            f.write('\n; Calibration frames\n')
            if 'bias' in scheduled_table['obj'].value:
                obj_bias = scheduled_table[scheduled_table['obj'] == 'bias']
                f.write("#COUNT {}\n#BIAS\n".format(obj_bias['count'][0])) # bias
                
            if 'dark' in scheduled_table['obj'].value:
                obj_dark = scheduled_table[scheduled_table['obj'] == 'dark']
                f.write("#COUNT {}\n#INTERVAL {}\n#DARK\n".format(obj_dark['count'][0], obj_dark['exptime'][0])) # dark
            
            if 'duskflat' in scheduled_table['obj'].value:
                f.write('#DUSKFLATS\n') # duskflat
            

            # Image frames
            f.write('\n; Start of the evening\n')
            f.write('#WAITUNTIL 1, {}\n'.format(time_astronomical_twilight_form))
            #f.write('#AFINTERVAL {}\n'.format(self.config['SCHEDULER_AFTIME'])) # AF every 3 hrs
            
            f.write('\n; Targeting\n')
            target_table = scheduled_table[scheduled_table['id'] != '']
            for target in target_table:
                if target['obj'] == 'autofocus':
                    f.write('\n#AUTOFOCUS\n')
                else:
                    name = target['obj']
                    ra = target['ra']
                    dec = target['dec']
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
            
            if 'dawnflat' in scheduled_table['obj'].value:
                f.write('#DAWNFLATS\n') # dawnflat
                
            if shutdown:
                f.write('\n; Closing\n')
                f.write('#SHUTDOWN\n')
            '''
            with open(os.path.dirname(self.rtspath)+'/autoflat_script.txt', 'w') as f_flat:
                f_flat.write('\n#WAITUNTIL 1, {}\n\n'.format(time_astronomical_dawn_form))
                f_flat.write('#DAWNFLATS\n')
                f_flat.write('\n; Closing\n')
                f_flat.write('#SHUTDOWN\n')            
            '''

    def write_ACPscript_RASA36(self,
                               filename_prefix : str = '',
                               savepath = './ACPscript/'
                               ):
        """
        Writes the observation plan in ACP script format for RASA36 telescope.

        Parameters
        ==========
        1. filename_prefix : str, optional, default=''
            A prefix to add to the filename when creating a new file.
        2. savepath : str, optional, default='./'
            The path to the directory where the file will be saved.

        Returns
        =======
        None
        """

        os.makedirs(savepath, exist_ok = True)
        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form
        
        # Set observing time 
        time_civil_twilight = self.civil_night[0]
        time_civil_dawn = self.civil_night[1]
        time_astronomical_twilight = self.tonight[0]
        time_astronomical_dawn = self.tonight[1]
        
        time_civil_twilight_form = dt_to_form(time_civil_twilight.datetime)
        time_civil_dawn_form = dt_to_form(time_civil_dawn.datetime)
        time_astronomical_twilight_form = dt_to_form(time_astronomical_twilight.datetime)
        time_astronomical_dawn_form = dt_to_form(time_astronomical_dawn.datetime)
        filename = f'{savepath}{filename_prefix}ACP_{(self.midnight.datetime).strftime("%y%m%d_%H%m")}_{self.name_telescope}.txt'
        
        with open(filename, 'w') as f:
            
            # Preparation
            f.write('; Prepare for the observation\n')
            f.write('#WAITUNTIL 1, {}\n'.format(time_civil_twilight_form))
            
            # Cooling
            if self.config['SCHEDULER_CCDCOOL']:
                f.write('\n; Cooling\n')
                f.write('#CHILL {}, {}\n'.format(self.config['SCHEDULER_CCDTEMP'], 1.0))
            
            # Calibration frames
            f.write('\n; Calibration frames\n')
            if 'bias' in scheduled_table['obj'].value:
                obj_bias = scheduled_table[scheduled_table['obj'] == 'bias']
                f.write("#COUNT {}\n#BIAS\n".format(obj_bias['count'][0])) # bias
                
            if 'dark' in scheduled_table['obj'].value:
                obj_dark = scheduled_table[scheduled_table['obj'] == 'dark']
                f.write("#COUNT {}\n#INTERVAL {}\n#DARK\n".format(obj_dark['count'][0], obj_dark['exptime'][0])) # dark
            
            if 'duskflat' in scheduled_table['obj'].value:
                f.write('#DUSKFLATS\n') # duskflat
            

            # Image frames
            f.write('\n; Start of the evening\n')
            f.write('#WAITUNTIL 1, {}\n'.format(time_astronomical_twilight_form))
            #f.write('#AFINTERVAL {}\n'.format(self.config['SCHEDULER_AFTIME'])) # AF every 3 hrs
            
            f.write('\n; Targeting\n')
            target_table = scheduled_table[scheduled_table['id'] != '']
            for target in target_table:
                if target['obj'] == 'autofocus':
                    f.write('\n#AUTOFOCUS\n')
                else:
                    name = target['obj']
                    ra = target['ra']
                    dec = target['dec']
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
            
            if 'dawnflat' in scheduled_table['obj'].value:
                f.write('#DAWNFLATS\n') # dawnflat

            '''
            with open(os.path.dirname(self.rtspath)+'/autoflat_script.txt', 'w') as f_flat:
                f_flat.write('\n#WAITUNTIL 1, {}\n\n'.format(time_astronomical_dawn_form))
                f_flat.write('#DAWNFLATS\n')
                f_flat.write('\n; Closing\n')
                f_flat.write('#SHUTDOWN\n')            
            '''
     
    def write_ACPscript_LSGT(self,
                             filename_prefix : str = '',
                             savepath = './ACPscript/',
                             period_script : int = 3 # hour
                             ):
        """
        Writes the observation plan in ACP script format for LSGT telescope.

        Parameters
        ==========
        1. filename_prefix : str, optional, default=''
            A prefix to add to the filename when creating a new file.
        2. savepath : str, optional, default='./'
            The path to the directory where the file will be saved.
        3. period_script : int, optional, default=3
            The number of hours in each period of the script.

        Returns
        =======
        None
        """
        
        os.makedirs(savepath, exist_ok = True)
        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form
        def ut_to_lt(utctime):
            return self.observer_tcspy.localtime(utctime)
        
        # Set observing time 
        time_astronomical_dawn = self.tonight[1]
        time_astronomical_dawn_form = dt_to_form(time_astronomical_dawn.datetime)
        
        target_table = scheduled_table[scheduled_table['id'] != '']
        
        time_start = target_table[0]['time_start']
        time_tot = (target_table[-1]['time_end'] - target_table[0]['time_start']).total_seconds()/3600
        time_grid = Time(time_start) + period_script * np.arange(time_tot//period_script) * u.hour
        
        for timecut in time_grid:
            if timecut != time_grid[-1]:
                idx = (timecut < target_table['time_start']) & ( target_table['time_start'] < timecut + period_script * u.hour)
                split_table = target_table[idx]
                filename = f'{savepath}{filename_prefix}ACP_{(ut_to_lt(split_table[0]["time_start"])).strftime("%y%m%d_%H%m")}_{period_script}hour_{self.name_telescope}.txt'
                with open(filename, 'w') as f:
                    f.write('; Prepare for the observation\n')
                    f.write('\n; Targeting\n')
                    for target in split_table:
                        name = target['obj']
                        ra = target['ra']
                        dec = target['dec']
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
            else:
                idx = (timecut < target_table['time_start'])
                split_table = target_table[idx]
                filename = f'{savepath}{filename_prefix}ACP_{(ut_to_lt(split_table[0]["time_start"])).strftime("%y%m%d_%H%m")}_{period_script}Hour_{self.name_telescope}.txt'
                with open(filename, 'w') as f:
                    f.write('; Prepare for the observation\n')
                    f.write('\n; Targeting\n')
                    for target in split_table:
                        name = target['obj']
                        ra = target['ra']
                        dec = target['dec']
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

    def write_NINAscript(self,
                         filename_prefix : str = '',
                         savepath = './'
                         ):
        pass
        
        
    def _match_targetDB_format(self,
                               default_exptime : int = 120,
                               default_binning : int = 1,
                               default_count : int = 5):
        # Match format of the targetlist
        if isinstance(self._target_db, pd.core.frame.DataFrame):
            data = Table().from_pandas(self._target_db)
        elif isinstance(self._target_db, Table):
            data= self._target_db.copy()
        elif isinstance(self._target_db, dict):
            header = self._target_db['header']
            values = self._target_db['value']
            data_pd = pd.DataFrame(values, columns = header)
            data = Table().from_pandas(data_pd)
        
        # Match the columns of the targetlist
        for column in ['obj', 'ra', 'dec']:
            if not column in data.keys():
                raise ValueError(f'"{column}" does not exist in the data')
        default_values = dict(filter = self.config['SCHEDULER_DEFAULTFILT'], exptime = self.config['SCHEDULER_DEFAULTEXP'], count = self.config['SCHEDULER_DEFAULTCOUNT'], binning = self.config['SCHEDULER_DEFAULTBIN'], priority = self.config['SCHEDULER_DEFAULTPRIOR'], note ='' )
        for column in ['filter', 'exptime', 'count', 'binning', 'priority', 'note']:
            if not column in data.keys():
                data[column] = default_values[column]
        try:
            data['priority'] = data['priority'].astype(float)
        except:
            data['priority'] = '1'
            data['priority'].astype(float)
        data['filter'] = data['filter'].astype('str')    
        data['exptime'] = data['exptime'].astype('str')
        data['count'] = data['count'].astype('str')
        data['binning'] = data['binning'].astype('str')
        for i in range(len(data)):
            exptimes = data[i]['exptime'].astype('str').split(',')
            counts = data[i]['count'].astype('str').split(',')
            binnings = data[i]['binning'].astype('str').split(',')
            filters = data[i]['filter'].astype('str').split(',')
            n_set = len(filters)
            if not len(exptimes) == n_set:
                exptime = exptimes[0]
                if not isinstance(exptime, float):
                    exptime = default_exptime
                exptimes = ','.join([str(exptime)] * n_set)
            else:
                exptimes = ','.join(exptimes)
            if not len(counts) == n_set:
                count = counts[0]
                if not isinstance(count, float):
                    count = default_count
                counts = ','.join([str(count)] * n_set)
            else:
                counts = ','.join(counts)
            if not len(binnings) == n_set:
                binning = binnings[0]
                if not isinstance(binning, float):
                    binning = default_binning
                binnings = ','.join([str(binning)] * n_set)
            else:
                binnings = ','.join(binnings)
            data[i]['exptime'] = exptimes
            data[i]['count'] = counts
            data[i]['binning'] = binnings   
        return data
    
    def _match_date_format(self):
        try:
            # Try to parse as an astropy time
            date = Time(self._date)
            return date
        except ValueError:
            pass

        # Try to parse as one of the date formats
        formats = ['%Y%m%d', '%Y/%m/%d', '%Y.%m.%d', '%y%m%d', '%y/%m/%d', '%y.%m.%d']
        for fmt in formats:
            try:
                date_dt = datetime.datetime.strptime(self._date, fmt)
                date = Time(date_dt)
                return date
            except ValueError:
                pass

        # If none of the above worked, raise an exception
        raise ValueError(f'Invalid date format: {self._date}')
    
    def _match_coord_format(self, ra, dec):
        if (isinstance(ra, float)) | (isinstance(ra, int)) | (isinstance(ra, str)):
            ra = str(ra) ; dec = str(dec)
            if (':' in ra) & (':' in dec):
                skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'icrs')
            elif ('h' in ra) & ('d' in dec):
                skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'icrs')
            elif (' ' in ra) & (' ' in dec):
                skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'icrs')
            else:
                skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg), frame = 'icrs')
        else:
            try:
                if (isinstance(ra[0], float)):
                    skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg), frame = 'icrs')
                else:
                    if (':' in ra[0]) & (':' in dec[0]):
                        skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'icrs')
                    elif ('h' in ra[0]) & ('d' in dec[0]):
                        skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'icrs')
                    elif (' ' in ra[0]) & (' ' in dec[0]):
                        skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'icrs')
            except:
                raise ValueError('Input format is not supported')
        return skycoord
    
    def _get_moon_radec(self,
                        date):
        return self.observer_tcspy.moon_radec(date)
    
    def _get_midnight(self):
        return Time((self.tonight[1].jd + self.tonight[0].jd)/2, format = 'jd')
        
    def _get_tonight(self,
                     horizon : int = -18):
        if self._entire_night:
            sunrisetime = self.observer_tcspy.sun_risetime(self.date, mode = 'next', horizon= horizon)
            sunsettime = self.observer_tcspy.sun_settime(utctimes = sunrisetime, mode = 'previous', horizon= horizon)
        else:
            sunrisetime = self.observer_tcspy.tonight(self.date, horizon= horizon)[1]
            sunsettime = self.observer_tcspy.tonight(self.date,  horizon= horizon)[0]
        return sunsettime, sunrisetime
    
    def _get_observable_target(self,
                               fraction_observable : float = 0.05):
        # fraction_observable means (observable time/all night), If 0.05, then the amount of observable time is roughly 30 minutes
        all_targets = self._get_target_astroplan(obs_tbl = self.all_target_tbl)
        observability_tbl = observability_table(constraints = self.constraints, observer = self.observer_astroplan, targets = all_targets, time_range = [self.tonight[0], self.tonight[1]], time_grid_resolution = 10 * u.minute)
        key = observability_tbl['fraction of time observable'] > fraction_observable
        obs_tbl = self.all_target_tbl[key]
        return obs_tbl
    
    def _get_target_tcspy(self,
                          obs_tbl : Table):
        all_targets = [mainTarget(name_telescope = self._name_telescope, observer= self.observer_tcspy, target_ra = target['ra'], target_dec= target['dec'], target_name= target['obj']) for target in obs_tbl]
        return all_targets

    def _get_target_astroplan(self, 
                              obs_tbl : Table):
        all_targets = [FixedTarget(self._match_coord_format(target['ra'], target['dec']), name = target['obj']) for target in obs_tbl]
        return all_targets

    def _update_tbl(self,
                    obs_tbl : Table):
        target_tcspy = self._get_target_tcspy(obs_tbl = obs_tbl)
        target_astroplan = self._get_target_astroplan(obs_tbl = obs_tbl)
        target_coords = self._match_coord_format(ra = obs_tbl['ra'], dec = obs_tbl['dec'])
        target_ids = [uuid.uuid4().hex for i in range(len(obs_tbl))]
        obs_tbl['target_tcspy'] = target_tcspy
        obs_tbl['target_astroplan'] = target_astroplan
        obs_tbl['ra_deg'] = target_coords.ra.value.round(5)
        obs_tbl['dec_deg'] = target_coords.dec.value.round(5)
        
        transittimelist = []
        risetimelist = []
        settimelist = []
        maxaltlist = []
        for target in target_tcspy:
            transittime = target.meridiantime(self.midnight, n_grid_points = 20).datetime
            risetime = target.risetime(transittime, n_grid_points = 20, mode = 'previous').datetime
            settime = target.settime(transittime, n_grid_points = 20, mode = 'next').datetime
            if transittime < self.tonight[0]:
                maxalt = target.altaz(self.tonight[0]).alt.value
            elif transittime > self.tonight[1]:
                maxalt = target.altaz(self.tonight[1]).alt.value
            else:
                maxalt = target.altaz(transittime).alt.value
            transittimelist.append(transittime)
            risetimelist.append(risetime)
            settimelist.append(settime)
            maxaltlist.append(maxalt)
            
        obs_tbl['rise'] = risetimelist
        obs_tbl['transit'] = transittimelist
        obs_tbl['set'] = settimelist
        obs_tbl['maxalt'] = maxaltlist
        obs_tbl['scheduled'] = False
        obs_tbl['id'] = target_ids
        return obs_tbl

    def _set_constraints(self,
                        TARGET_MINALT : float = None,
                        TARGET_MAXALT : float = None,
                        TARGET_MAX_SUNALT : float = None,
                        TARGET_MOONSEP : float = None,
                        TARGET_MAXAIRMASS : float = None,
                        **kwargs) -> list:
        constraint_all = []
        if (TARGET_MINALT != None) & (TARGET_MAXALT != None):
            constraint_altitude = AltitudeConstraint(min = TARGET_MINALT * u.deg, max = TARGET_MAXALT * u.deg, boolean_constraint = False)
            constraint_all.append(constraint_altitude)
        if (TARGET_MAX_SUNALT != None):
            constraint_atnight = AtNightConstraint(max_solar_altitude = TARGET_MAX_SUNALT * u.deg)
            constraint_all.append(constraint_atnight)
        if TARGET_MOONSEP != None:
            constraint_gallatitude = MoonSeparationConstraint(min = TARGET_MOONSEP * u.deg, max = None)
            constraint_all.append(constraint_gallatitude)
        if TARGET_MAXAIRMASS != None:
            constraint_gallatitude = AirmassConstraint(min = 1, max = TARGET_MAXAIRMASS, boolean_constraint = False)
            constraint_all.append(constraint_gallatitude)
        return constraint_all

    def _set_weight_score(self,
                          TARGET_WEIGHT_RELATIVE_ALT = 0.5,
                          TARGET_WEIGHT_ABSOLUTE_ALT = 0.5,
                          TARGET_WEIGHT_PRIORITY = 0.5,
                          TARGET_WEIGHT_SLEW = 0,
                          **kwargs):
        weight = dict()
        weight['TARGET_WEIGHT_RELATIVE_ALT'] = TARGET_WEIGHT_RELATIVE_ALT
        weight['TARGET_WEIGHT_ABSOLUTE_ALT'] = TARGET_WEIGHT_ABSOLUTE_ALT
        weight['TARGET_WEIGHT_PRIORITY'] = TARGET_WEIGHT_PRIORITY
        weight['TARGET_WEIGHT_SLEW'] = TARGET_WEIGHT_SLEW
        return weight
    
    def _calculate_tot_exptime(self, 
                               target):
        exptimes = target['exptime'].astype(str).split(',')
        counts = target['count'].astype(str).split(',')
        tot_exptime = 0
        for exptime, count in zip(exptimes, counts):
            tot_exptime += (float(exptime) * float(count))
        return tot_exptime

    def _calculate_tot_overhead(self,
                                target,
                                sepdist_slew : float = 60 #deg
                                ):
        filters = target['filter'].astype(str).split(',')
        counts = target['count'].astype(str).split(',')
        num_totimg = np.array(counts).astype(float).sum()
        num_filt = len(filters)
        time_readout = self.config['OVERHEAD_TIME_READOUT'] * num_totimg
        time_filtchange = self.config['OVERHEAD_TIME_FILTCHAN'] * num_filt
        time_slew = (sepdist_slew/self.config['OVERHEAD_SPEED_SLEWING']) * 3
        time_overhead = time_readout + time_filtchange + time_slew
        return time_overhead
    
    def _update_status_obs_tbl(self,
                               obs_tbl,
                               target):
        updated_tbl = obs_tbl.copy()
        idx = (target['id'] == updated_tbl['id'])
        updated_tbl['scheduled'][idx] = True
        return updated_tbl
    
    def _insert_schedule(self, 
                        start_time,
                        schedule_table,
                        target,
                        score = -999,
                        replace : bool = False):
        
        output_table = schedule_table.copy()
        start_time = Time(start_time)
        tot_exptime = self._calculate_tot_exptime(target = target)
        end_time = start_time + tot_exptime * u.s
        if replace:
            mid_time = output_table['time_start'] +datetime.timedelta(seconds = (schedule_table['time_end'][0] - schedule_table['time_start'][0]).total_seconds()/2)
            remove_idx = (mid_time < end_time) & (mid_time > start_time)
            output_table.remove_rows(remove_idx)
        all_values = dict()
        all_values['time_start'] = start_time.datetime
        all_values['time_end'] = end_time.datetime
        all_values['score'] = np.round(score, 2)
        for column in output_table.columns:
            if column in target.keys():
                all_values[column] = target[column]
        all_values_list = []
        for column in output_table.columns:
            if column in all_values.keys():
                all_values_list.append(all_values[column])
            else:
                all_values_list.append('0000')
        output_table.add_row(all_values_list)
        return output_table, end_time

    def _transitioner(self,
                      name : str,
                      exptime : str,
                      count : str,
                      filter : str = '',
                      id_ : str = None
                      ):
        
        if id_ == None:
            id_ = uuid.uuid4().hex
        transitioner = Table(names = ['obj', 'exptime', 'count', 'filter', 'id'], dtype = ['str', 'str', 'str', 'str', 'str'])
        transitioner.add_row([name, str(exptime), str(count), str(filter), str(id_)])
        return transitioner[0]
#%%
from astropy.io import ascii


#%% Usage
if __name__ == '__main__':
    from googlesheet import GoogleSheet
    from astropy.io import ascii
    #gs = GoogleSheet()
    #data = gs.get_sheet_data('ToO', format_ ='Table')
    # load table
    
    data_original = ascii.read('./test_IMSNG.dat', format = 'fixed_width')
    date = Time.now()
    data = '/Users/hhchoi1022/Desktop/SkyGridCatalog_KMTNet_90.csv'
    data_original = ascii.read(data)
    data_original.rename_column('id', 'obj')
    subsequent_days = 1
    file_prefix ='ToO_'
    for i in range(subsequent_days):
        obs = ObsScheduler(target_db = data_original, name_telescope = 'KCT', entire_night = False, duplicate_when_empty= False)
        obs.write_ACPscript_KCT(filename_prefix =file_prefix, shutdown = True) # for RASA36, shutdown = False
        obs.show(save= True, filename_prefix = file_prefix)
    
        
#%%
if __name__ == '__main__':
    data = Table()
    file_prefix ='IMSNG_'
    data['obj'] = data_original['obj']
    data['ra'] = data_original['ra']
    data['dec'] = data_original['dec']
    data['priority'] = data_original['prob_vol']
    date = Time('2023-04-22T00:36:23', format = 'isot')
    date = '20230422'
    #date = Time.now()
    #obs = ObsScheduler(target_db = data, name_telescope = 'CBNUO', entire_night = False, date = date)
    #obs.write_rts(scheduled_only= False, filename_prefix =file_prefix)    
    #obs.show(save= True, filename_prefix = file_prefix, savepath = './script/')
    obs = ObsScheduler(target_db = data, name_telescope = 'LSGT', entire_night = False, date = date, duplicate_when_empty= False)
    obs.write_ACPscript_LSGT(filename_prefix =file_prefix, period_script= 3) # for RASA36, shutdown = False
    obs.show(save= True, filename_prefix = file_prefix, savepath = './script/')
    obs = ObsScheduler(target_db = data, name_telescope = 'KCT', entire_night = False, date = date)
    obs.write_ACPscript_KCT(filename_prefix =file_prefix, shutdown = True) # for RASA36, shutdown = False
    obs.show(save= True, filename_prefix = file_prefix, savepath = './script/')
    obs = ObsScheduler(target_db = data, name_telescope = 'RASA36', entire_night = False, date = date)
    obs.write_ACPscript_RASA36(filename_prefix =file_prefix) # for RASA36, shutdown = False
    obs.show(save= True, filename_prefix = file_prefix, savepath = './script/')
        
# %%
