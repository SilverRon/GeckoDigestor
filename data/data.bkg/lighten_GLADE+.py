#%%
import numpy as np
from astropy.table import Table
import time
#%%
try:
	cat
except NameError:
	# cat = Table.read('../data/GLADE+.fits', format='fits', memmap=True)
	st = time.time()
	# cat = Table.read('../data/GLADE+.fits', memmap=True)
	cat = Table.read('../data/GLADE+.fits', format='fits')
	delt = time.time() - st
	print(f"Time to read GLADE+.fits: {delt:.3f} seconds")
# %%
n_total = len(cat)
print(f"Total number of objects: {n_total}")
# %%
distance_cut = 1000.
_cat = cat[cat['col33'] < distance_cut]
print(f"Number of objects with distance < {distance_cut:.1f} Mpc: {len(_cat)}/{len(cat)} ({1e2*len(_cat)/len(cat):1.1f}%)")
# %%
selected_columns = [
	'col1',
	'col8',
	'col9',
	'col10',
	'col29',
	'col30',
	'col32',
	'col33',
	'col34',
	'col36',
	'col37',
	'col38',
	'col39',
	'col40',
]
print(f"Number of selected columns: {len(selected_columns)}/{len(cat.colnames)}")
# %%
cat 