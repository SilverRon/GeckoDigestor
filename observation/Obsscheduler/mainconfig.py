#%%
# Written by Hyeonho Choi 2023.01
import glob, os
from astropy.io import ascii
from astropy.table import Table
import json
#%%
class mainConfig:
    def __init__(self,
                 name_telescope : str = 'KCT',
                 **kwargs):
        self.name_telescope = name_telescope
        configfilekey = os.path.dirname(__file__)+f'/config/{self.name_telescope}'+'/*.config',
        self._configfiles = glob.glob(configfilekey[0])
        if len(self._configfiles) == 0:
            print('No configuration file is found.\nTo make default configuration files, run tcspy.configuration.make_config')
        else:
            pass
            self.config = self._load_configuration()
            
    def _load_configuration(self):
        all_config = dict()
        for configfile in self._configfiles:
            with open(configfile, 'r') as f:
                config = json.load(f)
                all_config.update(config)
        return all_config

    def _initialize_config(self,
                           savedir = f'{os.path.dirname(__file__)}'):
        savepath = savedir + f'/config/{self.name_telescope}/'
        if not os.path.exists(savepath):
            os.makedirs(savepath, exist_ok= True)
        def make_configfile(dict_params, 
                            filename : str,
                            savepath = savepath):
            with open(savepath + filename, 'w') as f:
                json.dump(dict_params, f , indent = 4)
            print('New configuration file made : %s'%(savepath+filename))
        ###### ALL CONFIGURATION PARAMETERS(EDIT HERE!!!) ######
        observer_params = dict(OBSERVER_LONGITUDE = 127.97805, # Obstec = -70.7804,
                       OBSERVER_LATITUDE = 37.56583, # Obstec = -30.4704,
                       OBSERVER_ELEVATION = 140, # Obstec = 1580,
                       OBSERVER_TIMEZONE = 'Asia/Seoul', #'America/Santiago'
                       OBSERVER_NAME = 'Hyeonho Choi',
                       OBSERVER_OBSERVATORY = self.name_telescope
                       )
        target_params = dict(TARGET_MINALT = 30,
                             TARGET_MAXALT = 90,
                             TARGET_MAX_SUNALT = None,
                             TARGET_MOONSEP = 40,
                             TARGET_MAXAIRMASS = None,
                             TARGET_MERIDIAN_SEC = 1800,
                             TARGET_MAXSLEW = 100,
                             
                             TARGET_WEIGHT_RELATIVE_ALT = 0.5,
                             TARGET_WEIGHT_ABSOLUTE_ALT = 0,
                             TARGET_WEIGHT_PRIORITY = 0.5,
                             TARGET_WEIGHT_SLEW = 0.3
                             )
        overhead_params = dict(OVERHEAD_TIME_READOUT = 16,
                               OVERHEAD_TIME_FILTCHAN = 12,
                               OVERHEAD_TIME_AUTOFOCUS = 450,
                               OVERHEAD_SPEED_SLEWING = 5)
        
        scheduler_params = dict(SCHEDULER_AUTOFOCUS = True,
                                SCHEDULER_AFTIME = 3,
                                SCHEDULER_CCDCOOL = True,
                                SCHEDULER_CCDTEMP = -10,
                                
                                SCHEDULER_BIAS = True,
                                SCHEDULER_BIASCOUNT = 9,
                                
                                SCHEDULER_DARK = True,
                                SCHEDULER_DARKCOUNT = 9,
                                SCHEDULER_DARKEXPTIME = '120',
                                
                                SCHEDULER_FLAT_DUSK = False,
                                SCHEDULER_FLAT_DAWN = True,
                                SCHEDULER_FLATCOUNT = 9,
                                SCHEDULER_FLATEXPTIME = '120,120',
                                SCHEDULER_FLATFILTSET = 'g,r',
                                
                                SCHEDULER_DEFAULTFILT = 'g,r',
                                SCHEDULER_DEFAULTEXP = '120,120',
                                SCHEDULER_DEFAULTBIN = '1,1',
                                SCHEDULER_DEFAULTCOUNT = '5,5',
                                SCHEDULER_DEFAULTPRIOR = 1
                                )

        #make_configfile(observer_params, filename = f'Observer_{self.name_telescope}.config')
        #make_configfile(target_params, filename = f'Target_{self.name_telescope}.config')
        #make_configfile(overhead_params, filename = f'Overhead_{self.name_telescope}.config')
        make_configfile(scheduler_params, filename = f'Scheduler_{self.name_telescope}.config')

#%% Temporary running 
if __name__ == '__main__':
    telescope_list = ['KCT', 'RASA36', 'LSGT', 'CBNUO']
    for telescope in telescope_list:
        telescope = 'KCT'
        A = mainConfig(name_telescope = telescope)
        A._initialize_config()

#%%
# %%
