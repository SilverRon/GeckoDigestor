#%%
import os
from astropy.time import Time
from obsscheduler_V3 import ObsScheduler, ScriptMaker
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
#%%

dirlist = os.listdir('./Archive/')
for target in dirlist:
    name_project = 'GECKO'
    filename_prefix = 'GECKO_'
    date = Time.now()
    log_savepath = './log/'
    # ACP config
    ACP_savepath = f'./ACPscript/{target}/'
    # RTS config
    rts_savepath = f'./ACPscript/{target}/'
    n_target_for_each_timeslot = 2

    def get_isfile_and_data(target, name_telescope):


        host_file_key = f'./Archive/{target}/HostGalaxyCatalog_90.csv'
        is_host = os.path.isfile(host_file_key)
        grid_file_key = f'./Archive/{target}/SkyGridCatalog_{name_telescope}_90.csv'
        is_grid = os.path.isfile(grid_file_key)

        data_host = Table()
        data_targetted = Table()
        if is_host:
            data_host = ascii.read(host_file_key)
        if is_grid:
            data_targetted = ascii.read(grid_file_key)
        result = dict()
        result['host'] = dict()
        result['host']['exist'] = is_host
        result['host']['data'] = data_host
        result['grid'] = dict()
        result['grid']['exist'] = is_grid
        result['grid']['data'] = data_targetted
        return result
    
        
    
    ####################### 
    # KCT (ACP)
    #######################
    
    name_telescope = 'KCT'
    data = get_isfile_and_data(target = target, name_telescope= 'KCT')

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    if data['host']['exist']:
        scheduler_host = ObsScheduler(target_db= data['host']['data'],
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

    if data['grid']['exist']:
        scheduler_grid = ObsScheduler(target_db= data['grid']['data'],
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
    if data['grid']['exist']:
        plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
        plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    if data['host']['exist']:
        plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
        plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    
    
    #######################
    # RASA36 (ACP)
    #######################
    
    name_telescope = 'RASA36'
    data = get_isfile_and_data(target = target, name_telescope= name_telescope)

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    if data['host']['exist']:
        scheduler_host = ObsScheduler(target_db= data['host']['data'],
                                        date = date,
                                        name_project = name_project,
                                        name_telescope = name_telescope,
                                        entire_night = False)
        # Define target
        scriptmaker_host = ScriptMaker(scheduler_host)
        # Action
        scriptmaker_host.write_ACPscript_RASA36(filename_prefix= filename_prefix, savepath = ACP_savepath)
        scriptmaker_host.write_log(n_target = 300, sort_keyword = 'rank', filename_prefix= filename_prefix, savepath= ACP_savepath, format_ = 'ascii.fixed_width', return_ = False)
        scriptmaker_host.show(save = True, filename_prefix = filename_prefix, savepath = ACP_savepath)

    if data['grid']['exist']:
        scheduler_grid = ObsScheduler(target_db= data['grid']['data'],
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
    if data['grid']['exist']:
        plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
        plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    if data['host']['exist']:
        plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
        plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()

    
    #######################
    # LSGT (ACP)
    #######################
    name_telescope = 'LSGT'
    data = get_isfile_and_data(target = target, name_telescope= name_telescope)

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    if data['host']['exist']:
        scheduler_host = ObsScheduler(target_db= data['host']['data'],
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

    if data['grid']['exist']:
        scheduler_grid = ObsScheduler(target_db= data['grid']['data'],
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
    if data['grid']['exist']:
        plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
        plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    if data['host']['exist']:
        plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
        plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    
    
    ####################### 
    # LOAO (2x RTS)
    #######################
    
    name_telescope = 'LOAO'
    data = get_isfile_and_data(target = target, name_telescope= name_telescope)

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    if data['host']['exist']:
        scheduler_host = ObsScheduler(target_db= data['host']['data'],
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

    if data['grid']['exist']:
        scheduler_grid = ObsScheduler(target_db= data['grid']['data'],
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
    if data['grid']['exist']:
        plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
        plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    if data['host']['exist']:
        plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
        plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()


    ####################### 
    # CBNUO (2x RTS)
    #######################
    
    name_telescope = 'CBNUO'
    data = get_isfile_and_data(target = target, name_telescope= name_telescope)

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    if data['host']['exist']:
        scheduler_host = ObsScheduler(target_db= data['host']['data'],
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

    if data['grid']['exist']:
        scheduler_grid = ObsScheduler(target_db= data['grid']['data'],
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
    if data['grid']['exist']:
        plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
        plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    if data['host']['exist']:
        plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
        plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    
    ####################### 
    # SAO (2x RTS)
    #######################
    
    name_telescope = 'SAO'
    data = get_isfile_and_data(target = target, name_telescope= name_telescope)

    """
    Input
        1. targetted observation (HostGalaxyCatalog_90.csv)
        2. tiled obsevation (SkyGridCatalog_*_90.csv)
    Output 
        1. ACPscript
        2. log
    """
    if data['host']['exist']:
        scheduler_host = ObsScheduler(target_db= data['host']['data'],
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

    if data['grid']['exist']:
        scheduler_grid = ObsScheduler(target_db= data['grid']['data'],
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
    if data['grid']['exist']:
        plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
        plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    if data['host']['exist']:
        plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
        plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    
    
    #######################
    # CTIO (1x RTS)
    #######################
    
    name_telescope = 'KMTNet'
    data = get_isfile_and_data(target = target, name_telescope= name_telescope)

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
    if data['host']['exist']:
        scheduler_host = ObsScheduler(target_db= data['host']['data'],
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

    if data['grid']['exist']:
        scheduler_grid = ObsScheduler(target_db= data['grid']['data'],
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
    if data['grid']['exist']:
        plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
        plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    if data['host']['exist']:
        plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
        plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()


    ####################### 
    # SAAO (1x RTS)
    #######################
    
    name_telescope = 'KMTNet'
    data = get_isfile_and_data(target = target, name_telescope= name_telescope)

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
    if data['host']['exist']:
        scheduler_host = ObsScheduler(target_db= data['host']['data'],
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

    if data['grid']['exist']:
        scheduler_grid = ObsScheduler(target_db= data['grid']['data'],
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
    if data['grid']['exist']:
        plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
        plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    if data['host']['exist']:
        plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
        plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
    
    ######################## 
    # SSO (1x RTS)
    #######################
    
    name_telescope = 'KMTNet'
    data = get_isfile_and_data(target = target, name_telescope= name_telescope)

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
    if data['host']['exist']:
        scheduler_host = ObsScheduler(target_db= data['host']['data'],
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

    if data['grid']['exist']:
        scheduler_grid = ObsScheduler(target_db= data['grid']['data'],
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
    if data['grid']['exist']:
        plt.scatter(scheduler_grid.target.all['ra'], scheduler_grid.target.all['dec'], c = 'r', label ='grid', alpha = 0.3, s = 2)
        plt.scatter(scriptmaker_grid.result.scheduled['ra'], scriptmaker_grid.result.scheduled['dec'], c = 'b', label =f'scheduled[grid],{len(scriptmaker_grid.result.scheduled)}')
    if data['host']['exist']:
        plt.scatter(scheduler_host.target.all['ra'], scheduler_host.target.all['dec'], c = 'b', label ='host', alpha = 0.1, s = 2)
        plt.scatter(scriptmaker_host.result.scheduled['ra'], scriptmaker_host.result.scheduled['dec'], c = 'r', label =f'scheduled[host],{len(scriptmaker_host.result.scheduled)}')
    plt.legend()
