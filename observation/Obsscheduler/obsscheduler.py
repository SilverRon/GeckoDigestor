


#%%
import os
import sys
#dir_path = os.path.dirname(os.path.realpath(__file__))
#sys.path.append(os.path.join(dir_path, '..'))
from astropy.io import ascii
from mainobserver import mainObserver
from mainconfig import mainConfig
from maintarget import mainTarget
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astroplan import AltitudeConstraint, AirmassConstraint, MoonSeparationConstraint, AtNightConstraint
from astroplan import observability_table, FixedTarget
from astropy.time import Time
import datetime
import uuid
import matplotlib.pyplot as plt
import os
# %%

class ObsScheduler(mainConfig):
    """
    An observation scheduler class that generates observing schedules based on a target database, time constraints, and the telescope name.
    - History
    (23.04.22) Written by Hyeonho Choi 
    (23.04.23) Change in write_ACPscript_RASA36, write_ACPscript_KCT : change the order to apply #NOSOLVING for all targets
    (23.04.25) Scheduled table data format changed from Astropy.Table to Pandas.Dataframe
    (23.04.26) Added write_txt function
    (23.04.26) Change the output format of write_txt function
    (23.04.26) Change the structure of the configuration directory depending on the project (ToO, IMSNG)
    
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
    
    6. autofocus_at_init : bool (default=True)
        Whether to autofocus at the start of the plan

    Properties
    ==========
    1. target : The table for the information of all targets (including scheduled status)
    2. scheduls : The table for the scheduled plan

    Methods
    =======
    1. set_target(target_db : pd.Dataframe, Astropy.table.Table, Dict)
        Update the observation target list with new targets.
    2. set_date(date : Astropy.time.Time, datetime.datetime, str)
        Set the observation date for the scheduler.
    3. set_obs_info(name_telescope : str)
        Set up the observation information for the scheduler.
    1. scorer(obs_tbl: Table, telescope_altaz: SkyCoord, utctime: Time, datetime(Time.now()), duplicate: bool(False)) -> np.array
        Calculates scores for the given observation table based on various constraints and weights.
    2. scheduler() -> Table
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
                 project : str = 'ToO',
                 name_telescope = 'KCT',
                 entire_night : bool = False,
                 duplicate_when_empty : bool = False,
                 autofocus_at_init : bool = False
                 ):

        super().__init__(name_telescope = name_telescope, project = project)

        self._project = project
        self._name_telescope = name_telescope
        self._entire_night = entire_night
        self._duplicate = duplicate_when_empty
        self._autofocus_at_init = autofocus_at_init
        self.observer_tcspy = mainObserver(name_telescope= name_telescope, project = project)
        self.observer_astroplan = self.observer_tcspy._observer

        self.set_target(target_db = target_db)
        self.set_date(date = date)
        self.set_obs_info(name_telescope = name_telescope)
        self.scheduled_table = None
    
    @property
    def target(self):
        return self.target_all

    @property
    def schedule(self):
        return self.scheduled_table
    
    def set_target(self, 
                   target_db):
        """
        Update the observation target list with new targets.

        Parameters
        ----------
        target_db : pandas.DataFrame, astropy.table.Table, or dict
            A database of target objects, which must contain the columns 'obj', 'ra', and 'dec', as well as optional columns
            for 'filter', 'exptime', 'count', 'binning', 'weight', and 'note'.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If `target_db` is not a valid input format, or if a required column is missing.
            If `exptime`, `count`, and `binning` columns do not have the same length as the `filter` column.
            If the `priority` column is not float or int type.
        """
        
        # Match format of the targetlist
        target_db = target_db.copy()
        if isinstance(target_db, pd.core.frame.DataFrame):
            data = Table().from_pandas(target_db)
        elif isinstance(target_db, Table):
            data= target_db.copy()
        elif isinstance(target_db, dict):
            header = target_db['header']
            values = target_db['value']
            data_pd = pd.DataFrame(values, columns = header)
            data = Table().from_pandas(data_pd)
        else:
            raise ValueError("variable 'target_db' must be pd.Dataframe, astropy.table.Table, dict format")
        
        # Match the columns of the targetlist
        for column in ['obj', 'ra', 'dec']:
            if not column in data.keys():
                raise ValueError(f'"{column}" does not exist in the data')
        default_values = dict(filter = self.config['SCHEDULER_DEFAULTFILT'], 
                              exptime = self.config['SCHEDULER_DEFAULTEXP'], 
                              count = self.config['SCHEDULER_DEFAULTCOUNT'], 
                              binning = self.config['SCHEDULER_DEFAULTBIN'], 
                              priority = self.config['SCHEDULER_DEFAULTPRIOR'], 
                              note ='' )
        for column in ['filter', 'exptime', 'count', 'binning', 'weight', 'note']:
            if not column in data.keys():
                data[column] = default_values[column]
        
        
        
        
        
        try:
            data['weight'] = data['weight'].astype(float)
        except:
            raise ValueError("column 'weight' must be float or int")
        
        target_coords = self._match_coord_format(data['ra'], data['dec'])

        keys = data.keys()
        formatted_table = Table()

        formatted_table['obj'] = data['obj']
        formatted_table['ra_hms'] = [self._match_coord_to_string(coord)[0] for coord in target_coords]
        formatted_table['dec_dms'] = [self._match_coord_to_string(coord)[1] for coord in target_coords]
        formatted_table['filter'] = data['filter']
        formatted_table['exptime'] = data['exptime']
        formatted_table['count'] = data['count']
        formatted_table['binning'] = data['binning']
        formatted_table['weight'] = data['weight']
        formatted_table['note'] = data['note']
        formatted_table['scheduled'] = False
        for key in keys:
            formatted_table[key] = data[key]
        formatted_table['ra'] = np.round(target_coords.ra.value,4)
        formatted_table['dec'] = np.round(target_coords.dec.value,4)
        formatted_table['mainTarget'] = self._get_target_tcspy(obs_tbl = formatted_table)
        formatted_table['FixedTarget'] = self._get_target_astroplan(obs_tbl = formatted_table)
        formatted_table['id'] = [uuid.uuid4().hex for i in range(len(formatted_table))]

        for i in range(len(formatted_table)):
            exptimes = formatted_table[i]['exptime'].astype('str').split(',')
            counts = formatted_table[i]['count'].astype('str').split(',')
            binnings =formatted_table[i]['binning'].astype('str').split(',')
            filters = formatted_table[i]['filter'].astype('str').split(',')
            n_set = len(filters)
            if not len(exptimes) == n_set:
                raise ValueError("column 'exptime' must have same length as filters (comma separated input required) ")
            if not len(counts) == n_set:
                raise ValueError("column 'count' must have same length as filters (comma separated input required) ")
            if not len(binnings) == n_set:
                raise ValueError("column 'binning' must have same length as filters (comma separated input required) ")
        self.target_all = formatted_table.to_pandas()
    
    def set_date(self,
                 date):
        """
        Set the observation date for the scheduler.

        Parameters
        ----------
        date : str or astropy.time.Time
            The observation date, either as a string in a recognized date format, or as an Astropy Time object.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If `date` is not a valid input format.
        """
        
        if isinstance(date, str):
            formats = ['%y%m%d', '%y/%m/%d', '%y.%m.%d', '%Y%m%d', '%Y/%m/%d', '%Y.%m.%d']
            for fmt in formats:
                try:
                    date_dt = datetime.datetime.strptime(date, fmt)
                    date = Time(date_dt)
                    self.date = date
                    return 
                except: 
                    ValueError(f'Invalid date format: {date}')
        else:
            try:
                # Try to parse as an astropy time
                date = Time(date)
                self.date = date
                return 
            except: 
                ValueError(f'Invalid date format: {date}')

    def set_obs_info(self, name_telescope):
        """
        Set up the observation information for the scheduler.

        Parameters
        ----------
        name_telescope : str
            The name of the telescope used for observation.

        Returns
        -------
        None
        """
        
        super().__init__(name_telescope = name_telescope)
        self._name_telescope = name_telescope
        self.observer_tcspy = mainObserver(name_telescope= name_telescope)
        self.observer_astroplan = self.observer_tcspy._observer
        
        self.ast_allnight = self._get_allnight(horizon = -18)
        self.civil_allnight = self._get_allnight(horizon = 0)
        self.midnight = self._get_midnight()
        self.moon_radec = self._get_moon_radec(self.midnight)
        if self._entire_night:
            self.obs_night = self._get_allnight(horizon = -18)
        else:
            self.obs_night = self._get_obsnight(horizon = -18)
        self.target_observable = self._get_target_observable(0.05)
        if len(self.target_observable) == 0:
            raise RuntimeError('No observable target')
        target_df = self.target_observable.copy()
        risetime, transittime, settime, maxalt, maxaz = self._get_rts(obs_tbl = target_df)
        target_df['rise'] = risetime
        target_df['transit'] = transittime
        target_df['set'] = settime
        target_df['maxalt'] = maxalt
        ids = target_df['id']
        maintarget = target_df['mainTarget']
        fixedtarget = target_df['FixedTarget']
        target_df.drop(columns = ['mainTarget', 'FixedTarget', 'id'])
        target_df['FixedTarget'] = fixedtarget
        target_df['mainTarget'] = maintarget
        target_df['id'] = ids
        target_coords = self._match_coord_format(target_df['ra'].tolist(), target_df['dec'].tolist())
        all_target_altaz = self.observer_astroplan.altaz(self.midnight, target = target_coords)
        moon_altaz = self.observer_astroplan.altaz(self.midnight, target = self.moon_radec)
        moonsep = SkyCoord.separation(all_target_altaz, moon_altaz).value
        target_df['moonsep'] = moonsep
        
        self.target_observable = target_df
    
    def scheduler(self):
        # Set observing time 
        time_astronomical_twilight = self.ast_allnight[0]
        time_astronomical_dawn = self.ast_allnight[1]
        time_civil_twilight = self.civil_allnight[0]
        time_civil_dawn = self.civil_allnight[1]
        time_obs_start = self.obs_night[0]
        time_obs_end = self.obs_night[1]
        
        # Set transitioner
        target_bias = self._transitioner(name = 'bias', exptime = self.config['OVERHEAD_TIME_READOUT'], count = self.config['SCHEDULER_BIASCOUNT'], id_ ='')
        target_dark = self._transitioner(name = 'dark', exptime = self.config['SCHEDULER_DARKEXPTIME'], count = self.config['SCHEDULER_DARKCOUNT'], id_ ='')
        target_flat_dusk = self._transitioner(name = 'duskflat', exptime = self.config['SCHEDULER_FLATEXPTIME'], count = self.config['SCHEDULER_FLATCOUNT'], filter_ = self.config['SCHEDULER_FLATFILTSET'], id_ ='')
        target_flat_dawn = self._transitioner(name = 'dawnflat', exptime = self.config['SCHEDULER_FLATEXPTIME'], count = self.config['SCHEDULER_FLATCOUNT'], filter_ = self.config['SCHEDULER_FLATFILTSET'], id_ ='')
        
        obs_tbl = self.target_observable
        
        scheduled_table = self._initialize_schedule(tbl = obs_tbl)
        
        time_observation = time_obs_start
        # Preparation of the observation if observation start time is astronomical twilight
        if time_obs_start == time_astronomical_twilight:
            time_observation = time_civil_twilight
            # Dusk flat
            if self.config['SCHEDULER_FLAT_DUSK']:
                scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_flat_dusk, score = 1)
                target_overhead_flat = self._transitioner(name = f'overhead_flat', exptime = self._calculate_tot_overhead(target_flat_dusk), count = 1)
                scheduled_table,time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_overhead_flat)
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
        autofocus_time = time_observation
        if self._autofocus_at_init:
            autofocus_time = time_observation - 5 * u.hour
        telescope_altaz = obs_tbl.iloc[0].mainTarget.altaz(time_observation)
        target_best = None
        while time_observation < time_obs_end:
            # Autofocus
            if self.config['SCHEDULER_AUTOFOCUS']:
                if (time_observation - autofocus_time > self.config['SCHEDULER_AFTIME'] * u.hour) * (time_astronomical_dawn - time_observation > self.config['SCHEDULER_AFTIME'] * u.hour):
                    target_autofocus = self._transitioner(name = 'autofocus', exptime = self.config['OVERHEAD_TIME_AUTOFOCUS'], count = 1)
                    scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_autofocus, score = 1)
                    autofocus_time = time_observation
            # Target selection 
            target_best, telescope_altaz, target_score = self._scorer(obs_tbl, telescope_altaz, utctime = time_observation)
            
            if target_score == 0:
                if self._duplicate:
                    target_best, telescope_altaz, target_score = self._scorer(obs_tbl, telescope_altaz, utctime = time_observation, duplicate= True, order = 2)
                    if target_score == 0:
                        target_best = self._transitioner(name = f'empty_target', exptime = 600, count = 1)
                        target_score = -999
                else:
                    target_best = self._transitioner(name = f'empty_target', exptime = 600, count = 1)
                    target_score = -999
                
            scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_best, score = target_score)
            target_overhead = self._transitioner(name = f'overhead_{target_best["obj"]}', exptime = self._calculate_tot_overhead(target_best), count = 1)
            scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_overhead)
            obs_tbl = self._update_status_obs_tbl(obs_tbl = obs_tbl, target = target_best)
        
        if self.config['SCHEDULER_FLAT_DAWN']:
            scheduled_table, time_observation = self._insert_schedule(start_time = time_astronomical_dawn, schedule_table=scheduled_table, target = target_flat_dawn, score = 1)
            target_overhead_flat = self._transitioner(name = f'overhead_flat', exptime = self._calculate_tot_overhead(target_flat_dawn), count = 1)
            scheduled_table, time_observation = self._insert_schedule(start_time = time_observation, schedule_table=scheduled_table, target = target_overhead_flat)
        scheduled_table = scheduled_table[scheduled_table['score'] > 0 ]
        scheduled_table['scheduled'] = True
        self.target_all['scheduled'] = self.target_all.id.isin(scheduled_table.id).tolist()
        self.target_observable['scheduled'] = self.target_observable.id.isin(scheduled_table.id).tolist()
        self.scheduled_table = scheduled_table
        return scheduled_table
    
    def show(self,
             save : bool = False,
             filename_prefix : str = 'ToO_',
             savepath = './ACPscript/'
             ):
        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        target_table = scheduled_table.dropna(subset = ['ra'])
        target_table = Table().from_pandas(target_table)
        other_table = self.scheduled_table[self.scheduled_table['ra'].isna()]
        other_table = Table().from_pandas(other_table)
        n_obj = len(set(target_table['obj']))
        n_obs = len(set(target_table['id']))
        time_astronomical_dawn = self.ast_allnight[1]
        time_astronomical_twilight = self.ast_allnight[0]
        time_civil_dawn = self.civil_allnight[1]
        time_civil_twilight = self.civil_allnight[0]

        time_range_start, time_range_end = time_civil_twilight.datetime - datetime.timedelta(hours = 2), time_civil_dawn.datetime + datetime.timedelta(hours = 2)
        time_axis = np.arange(time_range_start, time_range_end, datetime.timedelta(minutes = 5))
        moon_altaz = self.observer_tcspy.moon_altaz(time_axis)
        sun_altaz = self.observer_tcspy.sun_altaz(time_axis)
        
        plt.figure(dpi = 300, figsize = (10, 4))
        plt.title(f'Observing plan[{self.name_telescope}]')
        for target in target_table:
            target_tcspy = target['mainTarget']
            target_altaz = target_tcspy.altaz(time_axis)
            obs_time_axis = np.arange(target['obs_start'].datetime, target['obs_end'].datetime, datetime.timedelta(minutes = 1))
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
        plt.text(time_civil_twilight.datetime-datetime.timedelta(minutes=140), 85, f'N_observation[N_target] = {(n_obs)}[{n_obj}]', fontsize = 10)
        plt.text(time_astronomical_twilight.datetime, 92, 'Ast.S.set', fontsize = 10)
        plt.text(time_civil_twilight.datetime-datetime.timedelta(minutes=00), 92, 'S.set', fontsize = 10)
        plt.text(time_civil_dawn.datetime-datetime.timedelta(minutes=00), 92, 'S.rise', fontsize = 10)
        plt.xlim(time_range_start - datetime.timedelta(hours = 1), time_range_end + datetime.timedelta(hours = 1))
        for target in other_table:
            plt.fill_betweenx([10,90], target['obs_start'].datetime, target['obs_end'].datetime, alpha = 0.1, color ='r')
        plt.fill_betweenx([-30,-31], target['obs_start'].datetime, target['obs_end'].datetime, alpha = 0.1, color ='r', label = 'Calib, Focus')
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
            
    def reset_scheduled_status(self):
        self.target_all['scheduled'] = False
        self.target_observable['scheduled'] = False
            
    def write_Gsheet(self,
                     spreadsheet_url : str,
                     authorize_jsonfile : str,
                     sheetname_prefix : str = '',
                     ):
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
    
    def write_txt(self,
                  scheduled_only : bool = False,
                  filename_prefix : str = 'ToO_',
                  savepath = './rts/',
                  format_ : str = 'ascii',
                  return_ : bool = True):
        def dt_to_form_hoursec(dt):
            try:
                form ='{:02}/{:02}:{:02}'.format(dt.day, dt.hour, dt.minute)
                return form    
            except:
                return None
            
        if scheduled_only:
            data_table = self.scheduled_table
        else:
            data_table = pd.merge(self.target_all, self.target_observable, how = 'left', on = None)
        data_table = data_table.dropna(subset = ['ra'])
        data_table = Table().from_pandas(data_table)
        risetime, transittime, settime, maxalt, maxaz = self._get_rts(obs_tbl = data_table, mode = 'tm')        
        target_coords = self._match_coord_format(data_table['ra'].tolist(), data_table['dec'].tolist())
        all_target_altaz = self.observer_astroplan.altaz(self.midnight, target = target_coords)
        moon_altaz = self.observer_astroplan.altaz(self.midnight, target = self.moon_radec)
        maxalt = np.array(maxalt).round(1)
        moonsep = np.array(SkyCoord.separation(all_target_altaz, moon_altaz).value).astype(int)
        data_table['maxalt'] = maxalt
        data_table['transit'] = transittime
        data_table['moonsep'] = moonsep
        
        filename = f'{savepath}{filename_prefix}{self._ut_to_lt(self.midnight.datetime).strftime("%y%m%d")}_{self.name_telescope}.log'
        # risetime_lt = []
        # transittime_lt = []
        # settime_lt = []
        
        # for target in data_table:
        #     try:
        #         rise_lt = dt_to_form_hoursec(self._ut_to_lt(target['rise'].datetime))
        #     except:
        #         rise_lt = None
        #     try:
        #         transit_lt = dt_to_form_hoursec(self._ut_to_lt(target['transit']))
        #     except:
        #         transit_lt = None
        #     try:
        #         set_lt = dt_to_form_hoursec(self._ut_to_lt(target['set'].datetime))
        #     except:
        #         set_lt = None
        #     risetime_lt.append(rise_lt)
        #     transittime_lt.append(transit_lt)
        #     settime_lt.append(set_lt)
        # data_table.remove_columns(['mainTarget', 'FixedTarget','id'])
        # data_table['rise'] = risetime_lt
        # data_table['transit'] = transittime_lt
        # data_table['set'] = settime_lt
        data_table.write(filename, format = format_, overwrite = True)
        if return_:
            return data_table
    
    def write_rts(self,
                  scheduled_only : bool = False,
                  filename_prefix : str = '',
                  savepath = './rts/'
                  ):
        os.makedirs(savepath, exist_ok = True)
        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()

        def dt_to_form_hoursec(dt):
            form ='{:02}:{:02}'.format(dt.hour, dt.minute)
            return form
        
        if scheduled_only:
            obs_info_pd = pd.merge(scheduled_table, self.target_observable, how = 'left', on = None)
        else:
            obs_info_pd = self.target_observable
        obs_info_pd = obs_info_pd.dropna(subset = ['ra'])
        obs_info = Table().from_pandas(obs_info_pd)
        obs_info.sort('ra')
        # Target
        all_target_radec = self._match_coord_format(ra = obs_info['ra'], dec = obs_info['dec'])
        all_target_altaz = self.observer_astroplan.altaz(self.midnight, target = all_target_radec)
        risetime_lt = [dt_to_form_hoursec(self._ut_to_lt(utctime.datetime)) for utctime in obs_info['rise']]
        transittime_lt = [dt_to_form_hoursec(self._ut_to_lt(utctime.datetime)) for utctime in obs_info['transit']]
        settime_lt = [dt_to_form_hoursec(self._ut_to_lt(utctime.datetime)) for utctime in obs_info['set']]
        moon_altaz = self.observer_astroplan.altaz(self.midnight, target = self.moon_radec)
        moonsep = SkyCoord.separation(all_target_altaz, moon_altaz).value.astype(int)
        # Moon
        moon_radec = self.observer_tcspy.moon_radec(self.midnight)
        moonphase = round(self.observer_tcspy.moon_phase(self.midnight),2)
        moon_ra_str = f'{int(moon_radec.ra.hms.h)}:{int(moon_radec.ra.hms.m)}:{round(moon_radec.ra.hms.s,2)}'
        moon_dec_str = f'{int(moon_radec.dec.dms.d)}:{np.abs(int(moon_radec.dec.dms.m))}:{round(np.abs(moon_radec.dec.dms.s),2)}'
        moon_set_sepration = int(self.config['TARGET_MOONSEP'])
        # Sun   
        sunrisetime = self._ut_to_lt(self.observer_tcspy.sun_risetime(self.date, mode = 'next', horizon= -18).datetime)
        sunsettime = self._ut_to_lt(self.observer_tcspy.sun_settime(utctimes = sunrisetime, mode = 'previous', horizon= -18).datetime)
        tonight = self._ut_to_lt(self.midnight.datetime).strftime('%Y-%m-%d')
        filename = f'{savepath}{filename_prefix}rts_{self._ut_to_lt(self.midnight.datetime).strftime("%y%m%d")}_{self.name_telescope}.txt'
        # Target table
        
        rts_table = Table()
        rts_table['name'] = obs_info['obj']
        rts_table['ra'] = obs_info['ra_hms']
        rts_table['dec'] = obs_info['dec_dms']
        rts_table['rise(LT)'] = risetime_lt
        rts_table['transit(LT)'] = transittime_lt
        rts_table['set(LT)'] = settime_lt
        rts_table['moon_dist(deg)'] = moonsep
        rts_table['weight'] = obs_info['weight']
        rts_table['note'] = obs_info['note']
        #rts_table['filter'] = obs_info['filter']
        #rts_table['exptime'] = obs_info['exptime']
        #rts_table['count'] = obs_info['count']
        #rts_table['binning'] = obs_info['binning']
        # Sorting
        #rts_table.sort('weight')
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
            f.write('='*110+'\n')
            for line in original_data:
                f.write(line)
            f.write('='*110+'\n')

    def write_ACPscript_KCT(self,
                            filename_prefix : str = 'ToO_',
                            savepath = './ACPscript/',
                            shutdown : bool = True,
                            ):
        os.makedirs(savepath, exist_ok = True)
        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form
        scheduled_table = Table().from_pandas(scheduled_table)
        # Set observing time 
        time_civil_twilight = self.civil_allnight[0] 
        time_civil_dawn = self.civil_allnight[1] 
        time_astronomical_twilight = self.obs_night[0]
        time_astronomical_dawn = self.obs_night[1]
        if self._entire_night:
            time_astronomical_twilight = self.ast_allnight[0]
            time_astronomical_dawn = self.ast_allnight[1]
        time_civil_twilight_form = dt_to_form(time_civil_twilight.datetime)
        time_civil_dawn_form = dt_to_form(time_civil_dawn.datetime)
        time_astronomical_twilight_form = dt_to_form(time_astronomical_twilight.datetime)
        time_astronomical_dawn_form = dt_to_form(time_astronomical_dawn.datetime)
        filename = f'{savepath}{filename_prefix}ACP_{(self.midnight.datetime).strftime("%y%m%d_%H%M")}_{self.name_telescope}.txt'
        
        with open(filename, 'w') as f:
            
            # Preparation
            f.write('; Prepare for the observation (KCT)\n')
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
            f.write('\n; Start of the observation\n')
            f.write('#WAITUNTIL 1, {}\n'.format(time_astronomical_twilight_form))
            #f.write('#AFINTERVAL {}\n'.format(self.config['SCHEDULER_AFTIME'])) # AF every 3 hrs
            
            f.write('\n; Targeting\n')
            target_table = scheduled_table[scheduled_table['id'] != '']
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
            
            if 'dawnflat' in scheduled_table['obj'].value:
                f.write('#DAWNFLATS\n') # dawnflat
                
            if shutdown:
                f.write('\n; Closing\n')
                f.write('#SHUTDOWN\n')

    def write_ACPscript_RASA36(self,
                               filename_prefix : str = 'ToO_',
                               savepath = './ACPscript/'
                               ):
        os.makedirs(savepath, exist_ok = True)
        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form
        scheduled_table = Table().from_pandas(scheduled_table)
        # Set observing time 
        time_civil_twilight = self.civil_allnight[0] 
        time_civil_dawn = self.civil_allnight[1] 
        time_astronomical_twilight = self.obs_night[0]
        time_astronomical_dawn = self.obs_night[1]
        if self._entire_night:
            time_astronomical_twilight = self.ast_allnight[0]
            time_astronomical_dawn = self.ast_allnight[1]
        time_civil_twilight_form = dt_to_form(time_civil_twilight.datetime)
        time_civil_dawn_form = dt_to_form(time_civil_dawn.datetime)
        time_astronomical_twilight_form = dt_to_form(time_astronomical_twilight.datetime)
        time_astronomical_dawn_form = dt_to_form(time_astronomical_dawn.datetime)
        filename = f'{savepath}{filename_prefix}ACP_{(self.midnight.datetime).strftime("%y%m%d_%H%M")}_{self.name_telescope}.txt'
        

        with open(filename, 'w') as f:
            
            # Preparation
            f.write('; Prepare for the observation (RASA36)\n')
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
            
            if 'dawnflat' in scheduled_table['obj'].value:
                f.write('#DAWNFLATS\n') # dawnflat

     
    def write_ACPscript_LSGT(self,
                             filename_prefix : str = 'ToO_',
                             savepath = './ACPscript/',
                             period_script : int = 3, # hour
                             ):
        os.makedirs(savepath, exist_ok = True)
        if self.scheduled_table is None:
            self.scheduler()
        scheduled_table = self.scheduled_table.copy()
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form

        scheduled_table = scheduled_table.dropna(subset = ['ra'])
        scheduled_table = Table().from_pandas(scheduled_table)
        
        # Set observing time 
        time_astronomical_dawn = self.obs_night[1]
        if self._entire_night:
            time_astronomical_dawn = self.ast_allnight[1]
        time_astronomical_dawn_form = dt_to_form(time_astronomical_dawn.datetime)
        
        target_table = scheduled_table
        
        time_start = target_table[0]['obs_start']
        time_tot = (target_table[-1]['obs_end'] - target_table[0]['obs_start']).value * 24
        time_grid = Time(time_start) + period_script * np.arange(time_tot//period_script) * u.hour
        
        for timecut in time_grid:
            if timecut != time_grid[-1]:
                idx = (timecut <= target_table['obs_start']) & ( target_table['obs_start'] < timecut + period_script * u.hour)
                split_table = target_table[idx]
                filename = f'{savepath}{filename_prefix}ACP_{(self._ut_to_lt(split_table[0]["obs_start"].datetime)).strftime("%y%m%d_%H%M")}_{period_script}hour_{self.name_telescope}.txt'
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
            else:
                idx = (timecut < target_table['obs_start'])
                split_table = target_table[idx]
                filename = f'{savepath}{filename_prefix}ACP_{(self._ut_to_lt(split_table[0]["obs_start"].datetime)).strftime("%y%m%d_%H%M")}_{period_script}Hour_{self.name_telescope}.txt'
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

    def write_NINAscript(self,
                         filename_prefix : str = '',
                         savepath = './'
                         ):
        pass

    def _scorer(self,
                obs_tbl : Table,
                telescope_altaz : SkyCoord,
                utctime : datetime or Time = None,
                duplicate : bool = False,
                order : int  = 1 # target with nth largest score will be returned 
                ):
        
        if utctime is None:
            utctime = Time.now()
        if not isinstance(utctime, Time):
            utctime = Time(utctime)
       
        all_target_radec = self._match_coord_format(ra = obs_tbl['ra'].tolist(), dec = obs_tbl['dec'].tolist())
        all_target_altaz = self.observer_astroplan.altaz(utctime, target = all_target_radec)
        all_target_alt = all_target_altaz.alt.value
        all_target_priority = obs_tbl['weight'].astype(float)
        all_target_hourangle_tmp = utctime.datetime - obs_tbl['transit']
        all_target_hourangle = np.array([target_hourangle.total_seconds() for target_hourangle in all_target_hourangle_tmp])
        target_sepdist = SkyCoord.separation(telescope_altaz, all_target_altaz).value
        
        score = np.ones(len(obs_tbl))
        # constraints
        constraint_moonsep = obs_tbl['moonsep'] > self.config['TARGET_MOONSEP']
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
        score_relative_alt = self.config['TARGET_WEIGHT_RELATIVE_ALT'] * all_target_alt / obs_tbl['maxalt']
        score_absolute_alt = self.config['TARGET_WEIGHT_ABSOLUTE_ALT'] * all_target_alt / 90
        score_priority = self.config['TARGET_WEIGHT_PRIORITY']* (all_target_priority / np.max(all_target_priority))
        if self._project == 'IMSNG':
            score_priority = self.config['TARGET_WEIGHT_PRIORITY']* (1- (all_target_priority / np.max(all_target_priority)))
        score_slewtime = self.config['TARGET_WEIGHT_SLEW'] * (target_sepdist / self.config['OVERHEAD_SPEED_SLEWING']) / np.max((target_sepdist / self.config['OVERHEAD_SPEED_SLEWING']))
        score_all = (score_relative_alt + score_absolute_alt + score_priority + score_slewtime) / 4
        score *= score_all
        if not np.sum(score) == 0:
            sorted_indices = np.argsort(np.array(score))
            idx_best = sorted_indices[-order]
            idx_best = np.argmax(np.array(score))
            target_best = obs_tbl.iloc[idx_best].to_frame().transpose()
            target_altaz = all_target_altaz[idx_best]
            target_score = score.iloc[idx_best]
            return target_best, target_altaz, target_score
        else:
            return None, telescope_altaz, 0

    def _ut_to_lt(self, utctime):
        return self.observer_tcspy.localtime(utctime)

    def _get_moon_radec(self,
                        date):
        return self.observer_tcspy.moon_radec(date)
    
    def _calculate_tot_exptime(self, 
                               target):
        
        exptimes = str(target['exptime'].tolist()[0]).split(',')
        counts = str(target['count'].tolist()[0]).split(',')
        tot_exptime = 0
        for exptime, count in zip(exptimes, counts):
            tot_exptime += (float(exptime) * float(count))
        return tot_exptime

    def _calculate_tot_overhead(self,
                                target,
                                sepdist_slew : float = 60 #deg
                                ):
        filters = str(target['filter'].tolist()[0]).split(',')
        counts = str(target['count'].tolist()[0]).split(',')
        num_totimg = np.array(counts).astype(float).sum()
        num_filt = len(filters)
        time_readout = self.config['OVERHEAD_TIME_READOUT'] * num_totimg
        time_filtchange = self.config['OVERHEAD_TIME_FILTCHAN'] * num_filt
        time_slew = (sepdist_slew/self.config['OVERHEAD_SPEED_SLEWING']) * 3
        time_overhead = time_readout + time_filtchange + time_slew
        return time_overhead
    
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
    
    def _match_coord_to_string(self, 
                               coord : SkyCoord):
        ra_hour = coord.ra.hms.h
        ra_minute = coord.ra.hms.m
        ra_seconds = coord.ra.hms.s
        dec_deg = coord.dec.dms.d
        dec_minute = np.abs(coord.dec.dms.m)
        dec_seconds = np.abs(coord.dec.dms.s)
        ra_string = '%02d:%02d:%02d'%(ra_hour, ra_minute, ra_seconds)
        dec_string = '%02d:%02d:%02d'%(dec_deg, dec_minute, dec_seconds)
        return ra_string, dec_string

    
    def _get_midnight(self):
        return Time((self.ast_allnight[1].jd + self.ast_allnight[0].jd)/2, format = 'jd')
        
    def _get_obsnight(self,
                      horizon : int = -18):
        sunrisetime = self.observer_tcspy.tonight(self.date, horizon= horizon)[1]
        sunsettime = self.observer_tcspy.tonight(self.date,  horizon= horizon)[0]
        return sunsettime, sunrisetime
    
    def _get_allnight(self,
                      horizon : int = -18):
        sunrisetime = self.observer_tcspy.sun_risetime(self.date, mode = 'next', horizon= horizon)
        sunsettime = self.observer_tcspy.sun_settime(utctimes = sunrisetime, mode = 'previous', horizon= horizon)
        return sunsettime, sunrisetime
    
    def _get_astroplan_constraints(self,
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
    
    def _get_target_observable(self,
                               fraction_observable : float = 0.05):
        # fraction_observable means (observable time/all night), If 0.05, then the amount of observable time is roughly 30 minutes
        constraints = self._get_astroplan_constraints(**self.config)
        alltargets = self.target_all['FixedTarget'].tolist()
        observability_tbl = observability_table(constraints = constraints, observer = self.observer_astroplan, targets = alltargets , time_range = [self.obs_night[0], self.obs_night[1]], time_grid_resolution = 10 * u.minute)
        key = observability_tbl['fraction of time observable'] > fraction_observable
        obs_tbl = self.target_all[key]
        return obs_tbl
    
    def _get_target_tcspy(self,
                          obs_tbl : Table):
        all_targets = [mainTarget(name_telescope = self._name_telescope, observer= self.observer_tcspy, target_ra = target['ra'], target_dec= target['dec'], target_name= target['obj']) for target in obs_tbl]
        return all_targets

    def _get_target_astroplan(self,
                              obs_tbl : Table):
        all_targets = [FixedTarget(self._match_coord_format(target['ra'], target['dec']), name = target['obj']) for target in obs_tbl]
        return all_targets

    def _get_rts(self,
                 obs_tbl,
                 mode : str  = 'rtsm'):
        transittimelist = []
        risetimelist = []
        settimelist = []
        maxaltlist = []
        maxazlist = []
        all_target = obs_tbl['mainTarget']
        for target in all_target:
            transittime = target.meridiantime(self.midnight, n_grid_points = 20)
            if 't' in mode:
                transittimelist.append(transittime.datetime)
            if 'r' in mode:
                risetime = target.risetime(transittime, n_grid_points = 20, mode = 'previous')
                risetimelist.append(risetime.datetime)
            if 's' in mode:
                settime = target.settime(transittime, n_grid_points = 20, mode = 'next')
                settimelist.append(settime.datetime)
            if 'm' in mode:
                if transittime < self.ast_allnight[0]:
                    maxaltaz = target.altaz(self.ast_allnight[0])
                    maxalt = np.round(maxaltaz.alt.value,2)
                    maxaz = np.round(maxaltaz.az.value,2)
                elif transittime > self.ast_allnight[1]:
                    maxaltaz = target.altaz(self.ast_allnight[1])
                    maxalt = np.round(maxaltaz.alt.value,2)
                    maxaz = np.round(maxaltaz.az.value,2)
                else:
                    maxalt = np.round(target.altaz(transittime).alt.value,2)
                    maxaz = np.round(target.altaz(transittime).az.value,2)
                maxaltlist.append(maxalt)
                maxazlist.append(maxaz)
        return risetimelist, transittimelist, settimelist, maxaltlist, maxazlist
    
    def _initialize_schedule(self,
                             tbl):
        empty_table = tbl.copy()
        empty_table = empty_table[0:0]
        return empty_table
    
    def _insert_schedule(self, 
                         start_time,
                         schedule_table,
                         target,
                         score = -999
                         ):
        output_table = schedule_table.copy()
        target_tmp = target.copy()
        start_time = Time(start_time)
        tot_exptime = self._calculate_tot_exptime(target = target)
        end_time = start_time + tot_exptime * u.s
        
        target_tmp['obs_start'] = start_time.datetime
        target_tmp['obs_end'] = end_time.datetime
        target_tmp['score'] = score
        if isinstance(target_tmp, pd.core.series.Series):
            target_df = target_tmp.to_frame().transpose()
        else:
            target_df = target_tmp
        output_table = pd.concat([output_table, target_df], axis = 0)
        return output_table, end_time
    
        
    def _update_status_obs_tbl(self,
                               obs_tbl,
                               target):
        updated_tbl = obs_tbl.copy()
        idx = (updated_tbl['id'] == target['id'].tolist()[0])
        updated_tbl.loc[idx, 'scheduled'] = True
        return updated_tbl
    
    def _transitioner(self,
                      name : str,
                      exptime : float or str,
                      count : float or str,
                      filter_ : str = '',
                      id_ : str = None
                      ):
        if id_ == None:
            id_ = uuid.uuid4().hex
        transitioner = pd.DataFrame(data = [[name, exptime, count, filter_, id_]], columns = ['obj', 'exptime', 'count', 'filter', 'id'])
        return transitioner
    
    

#%%

#%% Usage
if __name__ == '__main__':
    data = Table()
    file_prefix ='ToO'
    project = 'ToO'
    savepath = './ACPscript/'
    subsequent_day = 5
    
    data_original = ascii.read('/Users/hhchoi1022/Desktop/SkyGridCatalog_RASA36_90.csv')
    data_original.rename_column('prob_vol', 'weight')
    data_original.remove_column('id')
    #data_original = ascii.read('./test_IMSNG.dat', format = 'fixed_width')
    #data_original.rename_column('priority', 'weight')
    data['obj'] = data_original['obj']
    data['ra'] = data_original['ra']
    data['dec'] = data_original['dec']
    data['weight'] = data_original['weight']
    date = Time.now()
    
    date = Time.now()
    obs = ObsScheduler(target_db = data, name_telescope = 'CBNUO', entire_night = False, date = date)
    obs.write_rts(scheduled_only= False, filename_prefix =file_prefix)    
    obs.show(save= True, filename_prefix = file_prefix, savepath = './script/')
    obs = ObsScheduler(target_db = data, project = project, name_telescope = 'LSGT', entire_night = False, date = date, duplicate_when_empty= False)
    obs.write_rts(filename_prefix =file_prefix)
    obs.write_txt(scheduled_only = False, format_ = 'ascii.fixed_width', filename_prefix =file_prefix)
    obs.write_ACPscript_RASA36(filename_prefix =file_prefix) # for RASA36, shutdown = False   
    obs.write_ACPscript_KCT(filename_prefix =file_prefix, shutdown = True) # for RASA36, shutdown = False
    obs.write_ACPscript_LSGT(filename_prefix =file_prefix, period_script= 3, savepath = savepath) # for RASA36, shutdown = False
    obs.show(save= True, filename_prefix = file_prefix)
    #for i in range(subsequent_day):
    #date += 1* u.day
    obs = ObsScheduler(target_db = data, project = project, name_telescope = 'KCT', entire_night = True, date = date, duplicate_when_empty= True)
    obs.write_rts(filename_prefix =file_prefix)
    obs.write_txt(scheduled_only = False, format_ = 'ascii.fixed_width', filename_prefix =file_prefix)
    obs.write_ACPscript_KCT(filename_prefix =file_prefix, shutdown = True, savepath = './ACPscript/') # for RASA36, shutdown = False
    obs.show(save= True, filename_prefix = file_prefix)
    
    obs = ObsScheduler(target_db = data, project = project, name_telescope = 'RASA36', entire_night = True, date = date)
    obs.write_ACPscript_KCT(filename_prefix =file_prefix, shutdown = True) # for RASA36, shutdown = False
    obs.write_ACPscript_RASA36(filename_prefix =file_prefix, savepath =savepath) # for RASA36, shutdown = False
    obs.write_ACPscript_LSGT(filename_prefix =file_prefix, period_script= 3) # for RASA36, shutdown = False
    obs.write_rts(filename_prefix =file_prefix)
    obs.write_txt(scheduled_only = False, format_ = 'ascii.fixed_width', filename_prefix =file_prefix)
    obs.show(save= True, filename_prefix = file_prefix, savepath = savepath)

    obs = ObsScheduler(target_db = data, project = project, name_telescope = 'LOAO', entire_night = True, date = date)
    obs.write_ACPscript_KCT(filename_prefix =file_prefix, shutdown = True) # for RASA36, shutdown = False
    obs.write_ACPscript_RASA36(filename_prefix =file_prefix) # for RASA36, shutdown = False
    obs.write_ACPscript_LSGT(filename_prefix =file_prefix, period_script= 3) # for RASA36, shutdown = False
    obs.write_rts(filename_prefix =file_prefix)
    obs.write_txt(scheduled_only = False, format_ = 'ascii.fixed_width', filename_prefix =file_prefix)
    obs.show(save= True, filename_prefix = file_prefix, savepath = savepath)

    obs = ObsScheduler(target_db = data, project = project, name_telescope = 'CBNUO', entire_night = True, date = date)
    obs.write_ACPscript_KCT(filename_prefix =file_prefix, shutdown = True) # for RASA36, shutdown = False
    obs.write_ACPscript_RASA36(filename_prefix =file_prefix) # for RASA36, shutdown = False
    obs.write_ACPscript_LSGT(filename_prefix =file_prefix, period_script= 3) # for RASA36, shutdown = False
    obs.write_rts(filename_prefix =file_prefix)
    obs.write_txt(scheduled_only = False, format_ = 'ascii.fixed_width', filename_prefix =file_prefix)
    obs.show(save= True, filename_prefix = file_prefix, savepath = savepath)

    date = Time.now()
    obs = ObsScheduler(target_db = data, name_telescope = 'SAO', entire_night = True, date = date)
    obs.write_rts(scheduled_only= False, filename_prefix =file_prefix)    
    obs.write_rts(filename_prefix =file_prefix)
    obs.write_txt(scheduled_only = False, format_ = 'ascii.fixed_width', filename_prefix =file_prefix)
    obs.write_ACPscript_RASA36(filename_prefix =file_prefix) # for RASA36, shutdown = False   
    obs.write_ACPscript_KCT(filename_prefix =file_prefix, shutdown = True) # for RASA36, shutdown = False
    obs.write_ACPscript_LSGT(filename_prefix =file_prefix, period_script= 3, savepath = savepath) # for RASA36, shutdown = False
    obs.show(save= True, filename_prefix = file_prefix)
# %%
