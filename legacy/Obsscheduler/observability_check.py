#%%
from astropy.time import Time
from astropy.table import Table
from maintarget import mainTarget
import matplotlib.pyplot as plt
from mainobserver import mainObserver
import numpy as np
#%%
def check_observability(ra : str or float,
                        dec : str or float,
                        date  : Time = Time.now(),
                        check_telescopes : list = ['KCT', 'RASA36', 'LSGT', 'LOAO', 'SAO', 'CBNUO', 'KMTNet_CTIO', 'KMTNet_SSO', 'KMTNet_SAAO']
                        ):
    """
    Check observability of the target for given datetime & telescope

    Args:
        ra (str|float): target RA, in degree (supported format = 150.2233, '10:02:33')
        dec (str|float): target Dec, in degree (supported format = -33.2233, '-33:02:33')
        date (Time, optional): Datetime to check observability . Defaults to Time.now().
        check_telescopes (str|list, optional): Name of telescope(s). Defaults to ['KCT', 'RASA36', 'LSGT', 'LOAO', 'SAO', 'CBNUO', 'KMTNet_CTIO', 'KMTNet_SSO', 'KMTNet_SAAO'].
    """
    
    target = Table()
    target['obj'] = [objname]
    target['ra'] = ra
    target['dec'] = dec
    
    if isinstance(check_telescopes, str):
        observer = mainObserver(name_telescope = name_telescope, name_project = 'GECKO')
        maintarget = mainTarget(name_telescope = name_telescope, name_project= 'GECKO', observer = observer, target_ra = target['ra'], target_dec = target['dec'])
        maintarget.staralt(utctime = date)
        plt.title(name_telescope,f'[RA = {np.round(maintarget.ra,1)}, Dec = {np.round(maintarget.dec,1)}]')
    else:
        for name_telescope in check_telescopes:
            observer = mainObserver(name_telescope = name_telescope, name_project = 'GECKO')
            maintarget = mainTarget(name_telescope = name_telescope, name_project= 'GECKO', observer = observer, target_ra = target['ra'], target_dec = target['dec'])
            maintarget.staralt(utctime = date)
            plt.title(name_telescope+f'[RA = {np.round(maintarget.ra,1)}, Dec = {np.round(maintarget.dec,1)}]', fontsize = 10)
    
            

# %%
if __name__ == '__main__':
    objname = 'S230528ay'
    ra = 0.00
    dec = -30.00
    telescopes = ['KCT', 'RASA36', 'LSGT', 'LOAO', 'SAO', 'CBNUO', 'KMTNet_CTIO', 'KMTNet_SSO', 'KMTNet_SAAO']
    check_observability(ra = ra, dec = dec, check_telescopes= telescopes)
# %%
