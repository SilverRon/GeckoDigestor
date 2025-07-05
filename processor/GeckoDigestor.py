# %%
## Library
import os
# import glob
import time
import json
import astropy_healpix as ah
import numpy as np
from astropy.table import Table
from astropy.table import hstack
from astropy.coordinates import SkyCoord
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap.postprocess import crossmatch
from astropy import units as u
from astropy.io import ascii
import datetime
from astropy.time import Time
from base64 import b64decode
from io import BytesIO
import json
from pprint import pprint
from matplotlib.colors import LinearSegmentedColormap

# import tempfile
from astropy.table import Table
import astropy_healpix as ah
from astropy import units as u
from gcn_kafka import Consumer
import numpy as np
from astropy.io import fits
from ligo.skymap.io import write_sky_map
#	Location of the Sun
from astropy.coordinates import get_sun
from astropy.time import Time
import astropy.units as u
import healpy as hp
from healpy.newvisufunc import projview, newprojplot

import sys
sys.path.append('..')
# from config.config import *
# from src.AlertReceiver import *
from src.GeckoHelper import *
from util.util import *

#	plot setting
import matplotlib.pyplot as plt
import matplotlib as mpl
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "last_expr"

mpl.rcParams["axes.titlesize"] = 14
mpl.rcParams["axes.labelsize"] = 20
plt.rcParams['savefig.dpi'] = 500
plt.rc('font', family='serif')

import requests
def slack_bot(token, channel, text):
	response = requests.post("https://slack.com/api/chat.postMessage",
		headers={"Authorization": "Bearer "+token},
		data={"channel": channel,"text": text}
	)
	print(response)
# %%
import yaml

with open('../conf/config.yaml', 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)

# API 설정
api_config = config['API_CONFIG']
# GCN Kafka 설정
GCN_KAFKA_CONFIG = api_config['GCN_KAFKA']
client_id = GCN_KAFKA_CONFIG['client_id']
client_secret = GCN_KAFKA_CONFIG['client_secret']

# Slack 설정
SLACK_API_CONFIG = api_config['SLACK']
slack_token = SLACK_API_CONFIG['OAuth_Token']
slack_channel = SLACK_API_CONFIG['channel']

# GW 설정
gw_config = config['GW_SETTINGS']
TEMP_DIST_MEAN = gw_config['DISTANCE']['mean']
TEMP_DIST_STD = gw_config['DISTANCE']['std']

# 관측 설정
obs_config = config['OBSERVATION_SETTINGS']
# 경로 설정
path_out = obs_config['PATHS']['OUTPUT']
path_eventlog = obs_config['PATHS']['EVENT_LOG']
TEMPLATE_LC = obs_config['PATHS']['TEMPLATE_LC']

eventlogtbl = Table.read(path_eventlog, format='ascii.fixed_width')

# 데이터 경로
PATH_TO_GALAXY_CATALOG = obs_config['PATHS']['CATALOG']
PATH_TO_FACILITY = obs_config['PATHS']['FACILITY']

# 타일링 제한
NUMBER_GALAXY_LIMIT = obs_config['LIMITS']['GALAXIES']
NUMBER_TILING_LIMIT = obs_config['LIMITS']['TILES']

# 망원경 설정
telescopes4tile = obs_config['TELESCOPES']['enabled']
EXPTIME_LIMIT_SMNET = obs_config['LIMITS']['EXPOSURE']['SMNET']
EXPTIME_LIMIT_KMTNET = obs_config['LIMITS']['EXPOSURE']['KMTNET']

# 시간 설정
time_config = config['TIME_SETTINGS']
TIME_TO_SLEEP = time_config['SLEEP']

# 시각화 설정
vis_config = config['VISUALIZATION']
plot_config = vis_config['PLOT']
plt.rcParams['savefig.dpi'] = int(plot_config['DPI'])
mpl.rcParams["axes.titlesize"] = int(plot_config['TITLE_SIZE'])
mpl.rcParams["axes.labelsize"] = int(plot_config['LABEL_SIZE'])

# Debugging Mode
debug_mode = config['DEBUG_MODE']

cmap = LinearSegmentedColormap.from_list(
    'white_orchid_purple',
    [(0.0, 'w'),
     (0.5, 'orchid'),
     (1.0, 'purple')]
)

# slack
from slack_sdk import WebClient

# Slack 설정 로드
slack_config = config['API_CONFIG']['SLACK']
slack_token = slack_config['OAuth_Token']
slack_channel = slack_config['channel']

client = WebClient(token=slack_token)

# 이벤트 처리 설정
event_config = config['EVENT_PROCESSING']
DISTANCE_THRESHOLDS = event_config['DISTANCE']
PROBABILITY_THRESHOLDS = event_config['PROBABILITY_THRESHOLDS']
FARCUT = PROBABILITY_THRESHOLDS['FARCUT']

print("Configurations loaded:")
print(f"Path to galaxy catalog: {PATH_TO_GALAXY_CATALOG}")
print(f"Path to facility      : {PATH_TO_FACILITY}")
print(f"Number of galaxy limit: {NUMBER_GALAXY_LIMIT}")
print(f"Number of tiling limit: {NUMBER_TILING_LIMIT}")
print(f"Time to sleep         : {TIME_TO_SLEEP}")
print(f"Slack channel         : {slack_channel}")
print(f"Distance mean         : {TEMP_DIST_MEAN} Mpc")
print(f"Distance std          : {TEMP_DIST_STD} Mpc")
print(f"Telescopes enabled    : {telescopes4tile}")
print(f"SMNET exposure limit  : {EXPTIME_LIMIT_SMNET} sec")
print(f"KMTNET exposure limit : {EXPTIME_LIMIT_KMTNET} sec")
print(f"Plot DPI              : {plot_config['DPI']}")

print(f"Path to galaxy catalog: {PATH_TO_GALAXY_CATALOG}")
print(f"Path to facility      : {PATH_TO_FACILITY}")
print(f"Number of galaxy limit: {NUMBER_GALAXY_LIMIT}")
print(f"Number of tiling limit: {NUMBER_TILING_LIMIT}")
print(f"Time to sleep         : {TIME_TO_SLEEP}")

# %% [markdown]
# # 1. Receive the GW alert

# %% [markdown]
# - Online Test
# Attempt to reuse cached catalog in global namespace
try:
    cat  # check if 'cat' is already defined
except NameError:
    # Read the heavy FITS catalog only once
    cat = Table.read(PATH_TO_GALAXY_CATALOG, format='fits')

# Build SkyCoord from cached table
coordinates = SkyCoord(
    cat['col9'] * u.deg,
    cat['col10'] * u.deg,
    cat['col33'] * u.Mpc
)


print("=====             Start GeckoDigestor.py           =====\n")
print("===== GECKO prowls the cosmic savannah for GW prey =====\n")
print("""                                        
                       .. .             
                      .@%=@=.             |
                     .=*@@@%=           --+--
                     .+=#@@@%  .......    |
                        =@@= .*%@@@@@%*.
     .-=#+@@*.   .+@@@@@@@@@+*@@@@@@@@@@
    =*#@@@@@@#. -@@@@@@@@@@@@@@@@@@@@@@=
    .+@@@**@@@++@@@@@@@@@@@@@@@@@@@@%+. 
    .#--@+=@@@@@@@@@@@@@@@@@@- +%#-.    
       .. .%@@@@@@#=*%@%+@@%.           
     ..:+@@@@@@@@@@#-.  =@@#%..         
   .*@@@@@@@%+:#@@@@@=  -%@@@@#.        
 .#@@@#+-..     +=#@+     -@@@*         
:@@*.          .=#@@#+.                 
%@-            .*%%%-:.                 
@%                :*.                   
%@.                                     
:%@:.                                   
 .+@@@@%:                               
   .....                                
                                        """)

print(f"--------------------------------------------------------\n")
# %%

# 1) 전역 캐시 변수
_eventlog_mtime = 0
_eventlog_ids   = set()




def _refresh_eventlog():
    """
    Refresh the in-memory set of processed superevent_ids
    only if the log file has been modified since last read.
    """
    global _eventlog_mtime, _eventlog_ids

    try:
        mtime = os.path.getmtime(path_eventlog)
    except FileNotFoundError:
        # log file not yet created
        return

    # If the file is new or has changed, re-load
    if mtime != _eventlog_mtime:
        # Read the log table once
        tbl = Table.read(path_eventlog, format='ascii.fixed_width')
        # Cache the superevent_id column as a Python set
        _eventlog_ids = set(tbl['superevent_id'][tbl['processed']==False])
        _eventlog_mtime = mtime

def read_and_update_record(path_record, temp_dist_mean, temp_dist_std):
    """
    Read record.json, fill missing skymap-derived fields if needed,
    and determine most_probable_event.
    Returns (record_dict, record_updated_flag, most_probable_event_str).
    """
    # Read record.json
    with open(path_record, 'r') as f:
        record = json.load(f)

    # Fields we expect to fill from skymap
    required = ['ramax', 'decmax', 'area_90', 'distmean', 'diststd']
    missing = [f for f in required if f not in record]
    updated = False

    if missing:
        # Only decode & parse skymap if really needed
        skymap_str = record['event'].get('skymap')
        if skymap_str:
            # Decode and load into Table
            skymap = Table.read(BytesIO(b64decode(skymap_str)))

            # Compute most-probable pixel and its RA/Dec
            level, ipix = ah.uniq_to_level_ipix(
                skymap[np.argmax(skymap['PROBDENSITY'])]['UNIQ']
            )
            ra, dec = ah.healpix_to_lonlat(
                ipix, ah.level_to_nside(level), order='nested'
            )

            # Compute 90% credible area
            skymap.sort('PROBDENSITY', reverse=True)
            level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
            pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
            prob = pixel_area * skymap['PROBDENSITY']
            cumprob = np.cumsum(prob)
            idx90 = cumprob.searchsorted(0.9)
            area_90 = pixel_area[:idx90].sum().to(u.deg**2).value

            # Fill defaults for missing distance metadata
            distmean = skymap.meta.get('DISTMEAN', temp_dist_mean)
            diststd  = skymap.meta.get('DISTSTD',  temp_dist_std)

            # Update record dict
            record.update({
                'ramax':    ra.deg,
                'decmax':   dec.deg,
                'area_90':  area_90,
                'distmean': distmean,
                'diststd':  diststd
            })
            updated = True

    # Determine most probable event classification
    if record['event']['group'] == 'CBC':
        mpe = max(
            record['event']['classification'],
            key=record['event']['classification'].get
        )
    else:
        mpe = 'Burst'

    return record, updated, mpe

def load_skygrid(obs, base_dir='../data/skygrid'):
    """Load center+footprint tables and return a combined skygrid table."""
    center    = Table.read(f'{base_dir}/{obs}/displaycenter.txt',   format='csv')
    footprint = Table.read(f'{base_dir}/{obs}/displayfootprint.txt', format='csv')
    return hstack([center, footprint])


def process_tiles(skygrid_cat, skymap, simple_galcat,
                  confidence_limit, nside):
                #   confidence_limit, probkey, nside):
	"""
	1) Filter skygrid to tiles intersecting the 90% region.
	2) Add galactic coords, initialize n_hostgalaxy & probability.
	3) For each tile: count host galaxies & sum probabilities.
	4) Fill missing probabilities, sort, rank, compute confidence levels.
	Returns the processed Table.
	"""
	# 1) Select tiles in 90% region
	# idx = count_skymap_within_fov(skygrid_cat, skymap, confidence_limit)
	# idx = count_skymap_within_fov_fast(skygrid_cat, skymap, confidence_limit, nside)
	idx = count_skymap_within_fov_fast(skygrid_cat, skymap, confidence_limit)
	sel = skygrid_cat[idx]
	# if len(sel) == 0:
		# return sel

	# 2) Initialize columns
	c = SkyCoord(sel['ra'], sel['dec'], unit='deg', frame='icrs')
	sel['l'] = c.galactic.l
	sel['b'] = c.galactic.b
	sel['n_hostgalaxy'] = 0
	sel['n_healpix'] = 0
	sel["sum_p3d"] = 0.0
	sel["sum_p3d_x_stmass"] = 0.0

	# initialize new stats columns
	sel['mean_distmu']  = 0.0
	sel['median_distmu']= 0.0

	# Galaxy (x) Pixel Mapping
	theta_galaxy = np.radians(90 - simple_galcat['dec'].value)
	phi_galaxy   = np.radians(simple_galcat['ra'].value)
	galaxy_pix = hp.ang2pix(nside,
						theta_galaxy,
						phi_galaxy,
						nest=False)

	for i, row in enumerate(sel):
		# 1) 꼭짓점→벡터
		lon  = np.radians([row['ra1'],row['ra2'],row['ra3'],row['ra4']])
		lat  = np.radians([row['dec1'],row['dec2'],row['dec3'],row['dec4']])
		vecs = hp.ang2vec(0.5*np.pi - lat, lon)
		# 2) 폴리곤 픽셀 리스트
		tile_pix = hp.query_polygon(nside, vecs, inclusive=True)
		# 3) 교집합으로 은하 선택
		mask = np.isin(galaxy_pix, tile_pix)
		# if mask.sum() > 0:
			# print(i)
			# break
		sel['n_hostgalaxy'][i] = np.count_nonzero(mask)
		sel["sum_p3d_x_stmass"][i] = simple_galcat["prob_vol_x_stmass"][mask].sum()

	# Tile [x] Skymap
	# prepare sky map → pixel mapping once
	theta_skymap = np.radians(90.0 - skymap['DEC'])
	phi_skymap = np.radians(skymap['RA'])
	skymap_pix = hp.ang2pix(nside, theta_skymap, phi_skymap, nest=False)

	# 3) Compute per-tile stats
	for i, row in enumerate(sel):
		# tile corners → 3D vectors
		lon  = np.radians([row['ra1'], row['ra2'], row['ra3'], row['ra4']])
		lat  = np.radians([row['dec1'], row['dec2'], row['dec3'], row['dec4']])
		vecs = hp.ang2vec(0.5*np.pi - lat, lon)

		# get pixels inside this tile
		tile_pix = hp.query_polygon(nside, vecs, inclusive=True)

		# mask galaxies lying in those pixels
		m = np.isin(skymap_pix, tile_pix)
		sel['n_healpix'][i] = np.count_nonzero(m)

		# fill the three stats
		sel['sum_p3d'][i] = np.sum(skymap['PROBDENSITY'][m])
		sel['mean_distmu'][i] = np.mean(skymap['DISTMU'][m][np.isfinite(skymap['DISTMU'][m])])
		sel['median_distmu'][i] = np.median(skymap['DISTMU'][m][np.isfinite(skymap['DISTMU'][m])])

	# 4) Post-process: fill NaN, sort, rank, confidence 
	for probkey in ["sum_p3d", "sum_p3d_x_stmass"]:
		probs = sel[probkey][~np.isnan(sel[probkey])]
		sel[probkey] /= probs.sum()
		# valid = ~np.isnan(probs)
		# if np.any(valid):
			# sel[probkey][~valid] = np.min(probs[valid])
		# else:
			# sel[probkey] = 0.0

	#	Rank
	# sel['rank'] = np.arange(len(sel))
	# sel['rank_p3d_x_stmass'] = np.arange(len(sel))[np.argsort(sel["sum_p3d_x_stmass"])[::-1]]
	# sel.sort("sum_p3d", reverse=True)

	# 1) sum_p3d 내림차순 정렬
	order = np.argsort(sel['sum_p3d'])[::-1]
	sel = sel[order]  # 또는 sel.sort('sum_p3d', reverse=True)

	# 2) 순위 부여: 0이 가장 높은 확률
	sel['rank_p3d'] = np.arange(len(sel))

	order_x = np.argsort(sel['sum_p3d_x_stmass'])#[::-1]
	# 만약 테이블을 재정렬하지 않고, 별도 컬럼으로 등수만 기록하고 싶다면:
	ranks_x = np.empty(len(sel), dtype=int)
	ranks_x[order_x] = np.arange(len(sel))
	sel['rank_p3d_x_stmass'] = ranks_x

	cum   = np.cumsum(sel["sum_p3d"])
	total = cum[-1]
	sel['confidence'] = 0.0
	for level in [1.0, 0.99, 0.95, 0.9, 0.5]:
		sel['confidence'][cum <= total * level] = level

	return sel


def save_skygrid_catalog(sel, obs, record, probkey, confidence_limit, path_output):
    """
    Create table columns ('obj','note','weight','meta') and save as CSV.
    """
    name = os.path.join(path_output, f"SkyGridCatalog_{obs}_{record['superevent_id']}_{confidence_limit:.0%}.csv")
    sel['obj']    = record['superevent_id']
    sel['note']   = sel['#id']
    sel['weight'] = sel[probkey]
    sel.meta = {
        'obs': obs,
        'superevent_id': record['superevent_id'],
        'alert_type': record['alert_type'],
        'confidence': confidence_limit,
        'ordering': probkey
    }
    sel.write(name, format='csv', overwrite=True)


# def plot_tiling_map(sel, skymap, obs, record, path_output):
# 	"""
# 	Draw Healpix Mollweide map with tile edges overlaid.
# 	"""
# 	fig = plt.figure(figsize=(10, 4))
# 	hp.mollview(
# 		skymap,
# 		title=f"{obs} Tiling for {record['superevent_id']}",
# 		cbar=True, cmap=cmap, flip='astro',
# 		# coord=['E'], 
# 		unit=r'$P_{3D}$',
# 	)
# 	hp.graticule(dpar=30, dmer=30, coord='E', c='silver', alpha=0.75)

# 	# overlay tile edges
# 	ras_list, decs_list = DrawTiles(sel)
# 	# for ras, decs in zip(ras_list, decs_list):
# 	# 	hp.projplot(ras, decs, lonlat=True, c='g', lw=1, alpha=0.7)

# 	for ii, (ras, decs) in enumerate(zip(ras_list, decs_list)):
# 		if ii < 100: # Top 100
# 			hp.projplot(ras, decs, lonlat=True, c='lime', lw=1, alpha=1.0, zorder=999)
# 		else:
# 			hp.projplot(ras, decs, lonlat=True, c='k', lw=1, alpha=0.75, zorder=998)


# 	out = os.path.join(path_output, f"tiling_{obs}.png")
# 	plt.tight_layout()
# 	plt.savefig(out, dpi=100)
# 	plt.close(fig)

def plot_tiling_map(sel, skymap, obs, record, path_output):
	"""
	Draw Healpix Mollweide map with tile edges overlaid using projview.
	"""


	# projview 에 준 rot 값
	rot_lon, rot_lat = 180, 0  

	# 1) Figure 생성
	fig = plt.figure(figsize=(10, 4))

	# 2) projview로 맵 그리기
	projview(
		skymap,
		coord=['E'],
		graticule=True,
		graticule_labels=True,
		latitude_grid_spacing=30,
		longitude_grid_spacing=30,
		projection_type='mollweide',
		title=f"{obs} Tiling for {record['superevent_id']}",
		cmap=cmap,
		# cbar=True,
		flip='astro',
		unit=r'$P_{3D}$',
		# rot=(rot_lon, rot_lat),
	)

	# # 3) 타일 모서리 불러와서 newprojplot으로 오버레이
	ras_list, decs_list = DrawTiles(sel)
	for ii, (ras, decs) in enumerate(zip(ras_list, decs_list)):
		color = 'lime' if ii < 100 else 'k'
		alpha = 1.0 if ii < 100 else 0.75
		zord  = 999    if ii < 100 else 998
		newprojplot(
			ras, decs,
			lonlat=True,       # deg 단위 lon/lat
			linestyle='-',     # 실선
			linewidth=1,
			color=color,
			alpha=alpha,
			zorder=zord
		)


	# 3) 타일 모서리 불러와서 newprojplot으로 오버레이
	# ras_list, decs_list = DrawTiles(sel)
	# for ii, (ras, decs) in enumerate(zip(ras_list, decs_list)):
	# 	color = 'lime' if ii < 100 else 'k'
	# 	alpha = 1.0   if ii < 100 else 0.75
	# 	zord  = 999   if ii < 100 else 998

	# 	# NumPy array 로 변환
	# 	ras_arr  = np.array(ras)
	# 	decs_arr = np.array(decs)

	# 	# 회전 적용
	# 	ras_rot  = (ras_arr + rot_lon) % 360
	# 	decs_rot = decs_arr + rot_lat

	# 	newprojplot(
	# 		ras_rot, decs_rot,
	# 		lonlat=True,
	# 		linestyle='-',
	# 		linewidth=1,
	# 		color=color,
	# 		alpha=alpha,
	# 		zorder=zord
	# 	)

	# 4) 저장 및 종료
	out = os.path.join(path_output, f"tiling_{obs}.png")
	plt.tight_layout()
	plt.savefig(out, dpi=100)
	plt.close(fig)


def plot_cumulative_dist(sel, obs, probkey, path_output):
	"""
	Draw cumulative P3D plot for the selected tiles.
	"""
	fig = plt.figure(figsize=(6,4))
	cum = np.cumsum(sel[probkey]) / np.sum(sel[probkey])
	plt.plot(cum, '-', mfc='w', mew=2, ms=6, c='k')

	confidences = [0.99,0.95,0.9,0.5]
	linestyles = ['-','-.',':','--']
	colors = ['r', 'firebrick','indianred','lightcoral',]

	for lvl, ls, color in zip(confidences, linestyles, colors):
		n_tiles = np.interp(lvl, cum, np.arange(len(cum)))
		label = f"{int(lvl*100)}% ({int(n_tiles)})"
		plt.axhline(y=lvl, ls=ls, lw=2, alpha=0.7, label=label, c=color)

	plt.xlabel("Number of Tiles")
	plt.ylabel(r"Cumulative $\rm P_{3D}$")
	plt.legend(loc='lower right')
	plt.grid(ls='--', alpha=0.5)

	out = os.path.join(path_output, f"tiling_{obs}_cumulative_dist.png")
	plt.title(obs)
	plt.tight_layout()
	plt.savefig(out, dpi=100)
	plt.close(fig)

def plot_tiling_map_zoom(sel, skymap, obs, record, path_output,
                         ra_center, dec_center, fov_deg=10, xsize=300):
	"""
	Draw a zoomed-in gnomonic projection around (ra_center, dec_center).
	- sel          : selected skygrid Table (with tile edges)
	- skymap       : healpix probability map array
	- obs          : observatory name (for title/filename)
	- record       : record dict (for superevent_id)
	- path_output  : output directory
	- ra_center    : center RA in degrees
	- dec_center   : center Dec in degrees
	- fov_deg      : field-of-view in degrees
	- xsize        : width of the output image in pixels
	"""
	# compute resolution: degrees per pixel
	reso = fov_deg / xsize

	fig = plt.figure(figsize=(6,6))

	# gnomonic projection around the given center
	hp.gnomview(
		skymap,
		rot=(ra_center, dec_center, 0),  # center lon, lat, roll
		xsize=xsize,
		reso=reso,
		title=f"{obs} Zoom @({ra_center:.1f},{dec_center:.1f})",
		# cmap='CMRmap_r',
		cmap=cmap,
		# cbar=True,
		unit=r'$P_{3D}$',
		flip='astro',
	)
	# draw graticules at 2° interval (or adjust as you like)
	hp.graticule(dpar=15, dmer=15, coord='E', c='silver', alpha=0.75)

	# overlay tile edges (they get projected by projplot)
	ras_list, decs_list = DrawTiles(sel)
	for ii, (ras, decs) in enumerate(zip(ras_list, decs_list)):
		if ii < 100: # Top 100
			hp.projplot(ras, decs, lonlat=True, c='lime', lw=1, alpha=1.0, zorder=999)
		else:
			hp.projplot(ras, decs, lonlat=True, c='k', lw=1, alpha=0.75, zorder=998)

	# save and close    
	out = os.path.join(path_output,
						f"tiling_{obs}_zoom_{ra_center:g}_{dec_center:g}.png")
	plt.savefig(out, dpi=100, bbox_inches='tight', pad_inches=0.1)
	plt.close(fig)


def get_prob(skymap, ra, dec, max_level=29):
	max_nside = ah.level_to_nside(max_level)
	level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
	index = ipix * (2**(max_level - level))**2

	sorter = np.argsort(index)
	match_ipix = ah.lonlat_to_healpix(ra, dec, max_nside, order='nested')
	idx = sorter[np.searchsorted(index, match_ipix, side='right', sorter=sorter) - 1]
	p3d = skymap[idx]['PROBDENSITY'] * (np.pi / 180)**2
	return p3d


import numpy as np
import healpy as hp
import astropy_healpix as ah

def select_tiles_and_sum_prob(skygrid_cat, skymap_tbl, confidence_limit, nside_high):
    """
    1) multiresolution skymap의 UNIQ → (level, ipix_low)
    2) 각 low‐res 픽셀을 nside_high 해상도로 업샘플링
    3) 90% region 에 속하는 high‐res 픽셀 집합 ipix90_high 구성
    4) 각 타일 내부 픽셀 pix_tile 을 뽑아서
       (a) ipix90_high 과 교집합 여부로 타일 선택
       (b) 그 교집합 픽셀들의 PROBDENSITY 합산
    Returns:
      sel_idx: 선택된 타일의 인덱스 리스트
      sum_probs: 그 타일별 누적 PROBDENSITY 값 리스트
    """
    # 1) UNIQ → (level, ipix_low)
    level, ipix_low = ah.uniq_to_level_ipix(skymap_tbl['UNIQ'])
    nside_low = ah.level_to_nside(level)
    # 2) 업샘플 비율
    scale = nside_high // nside_low
    #  low‐pix 하나당 대응되는 high‐pix 첫 인덱스
    ipix_high_base = ipix_low * (scale**2)

    # PROBDENSITY 배열
    prob_density = skymap_tbl['PROBDENSITY']

    # 3) 90% 영역 high‐res 픽셀 집합
    mask90 = skymap_tbl['CUMPROBDENSITY'] < confidence_limit
    ipix90_high = set(np.unique(ipix_high_base[mask90]))

    sel_idx   = []
    sum_probs = []

    for i, row in enumerate(skygrid_cat):
        # 타일 꼭짓점 (rad)
        lons = np.radians([row['ra1'], row['ra2'], row['ra3'], row['ra4']])
        lats = np.radians([row['dec1'],row['dec2'],row['dec3'],row['dec4']])
        vecs = hp.ang2vec(0.5*np.pi - lats, lons)

        # 4-1) nside_high 해상도로 타일 내부 픽셀 인덱스
        pix_tile = hp.query_polygon(
            nside_high, vecs,
            inclusive=True, nest=True
        )

        # 4-2) 90% region 과 교집합이 있으면 이 타일을 채택
        if ipix90_high.intersection(pix_tile):
            sel_idx.append(i)

            # 4-3) 교집합 픽셀에 대응하는 PROBDENSITY 합산
            #   ipix_high_base 에 대응하는 픽셀들이 pix_tile 에 속하는지 검사
            mask_tile = np.isin(ipix_high_base, pix_tile)
            sum_probs.append(prob_density[mask_tile].sum())

    return np.array(sel_idx, dtype=int), np.array(sum_probs)

def clean_skygrid(skygrid_cat):
    """
    skygrid_cat 의 ra1..ra4, dec1..dec4 컬럼 중
    하나라도 masked(혹은 nan)인 행은 통째로 제거합니다.
    """
    n = len(skygrid_cat)
    bad = np.zeros(n, dtype=bool)

    # astropy MaskedColumn 인지, numpy array 인지 체크하며
    # mask 속성 혹은 np.isnan 로 bad 를 표시
    for col in ['ra1','ra2','ra3','ra4','dec1','dec2','dec3','dec4']:
        coldata = skygrid_cat[col]
        if hasattr(coldata, 'mask'):
            bad |= coldata.mask
        else:
            bad |= np.isnan(coldata)

    # bad==True 인 행을 모두 제거
    clean_tbl = skygrid_cat[~bad]
    n_removed = bad.sum()
    print(f"→ skygrid_cat 에서 {n_removed} 개의 결측 행 제거됨.")
    return clean_tbl


# %%
start_time = time.time()
while True:	
	# --- elapsed time 출력 (overwrite) ---
	elapsed = datetime.timedelta(seconds=int(time.time() - start_time))
	print(f"\rElapsed: [{elapsed.days:02d}:{elapsed.seconds//3600:02d}:"
			f"{(elapsed.seconds//60)%60:02d}:{elapsed.seconds%60:02d}]", end='', flush=True)

	# --- 1) 로그 갱신 ---
	# print(f"Refreshing event log...")
	# _refresh_eventlog()

	# # --- 2) 디렉터리 스캔 (os.scandir이 glob보다 빠릅니다) ---
	# # print(f"Scanning directory...")
	# eventlist = []
	# with os.scandir(path_out) as it:
	# 	for entry in it:
	# 		# only include directories, skip RETRACTION
	# 		if entry.is_dir() and 'RETRACTION' not in entry.name:
	# 			eventlist.append(entry.name)

	# # --- 3) 새로운 이벤트만 필터링 (set 차집합) ---
	# #    membership test in set is O(1)
	# new_events = [e for e in eventlist if e not in _eventlog_ids]

	eventlist = []
	with os.scandir(path_out) as it:
		for entry in it:
			# only include directories, skip RETRACTION
			if entry.is_dir() and 'RETRACTION' not in entry.name:
				eventlist.append(entry.name)

	# --- 3) 새로운 이벤트만 필터링 (set 차집합) ---
	#    membership test in set is O(1)
	new_events = []
	# for e in eventlist:
	# 	# Extract superevent_id and alert_type from directory name
	# 	parts = e.split('_')
	# 	superevent_id = '_'.join(parts[:-1])
	# 	alert_type = parts[-1]
		
	# 	# Check if event is in event log and not processed
	# 	if superevent_id in _eventlog_ids:
	# 		row = eventlogtbl[_eventlogtbl['superevent_id'] == superevent_id]
	# 		if len(row) > 0 and not row['processed'][0]:
	# 			new_events.append(e)
	
	eventlogtbl = Table.read(f"{path_out}/event.log", format='ascii.fixed_width')
	for e in eventlist:
		indx_not_processed = eventlogtbl['processed'] == False
		not_processed_events = eventlogtbl['output_dir'][indx_not_processed]
		if e in not_processed_events:
			new_events.append(e)

	if new_events:
		print(f"\nNew {len(new_events)} events found:", new_events)
		# TODO: new_events 처리 로직 호출
		# 예) for ev in new_events: handle_event(ev)

	# sleep or other logic
	time.sleep(TIME_TO_SLEEP)

	if len(new_events):
	# for ee, event in enumerate(new_events):
		event = new_events[0]
		st = time.time()
		print(f"{event} DETECTED")
 
		#	Path
		path_output = f"{path_out}/{event}"
		path_record = f"{path_output}/record.json"

		#	Read the record
		if not os.path.exists(path_record):
			continue  # skip if no record.json

		#	Define the path to save the results

		print(f"Processing {path_record}...")
		record, updated, most_probable_event = read_and_update_record(
			path_record, TEMP_DIST_MEAN, TEMP_DIST_STD
		)

		skymap_str = record['event']['skymap']
		skymap = Table.read(BytesIO(b64decode(skymap_str)))

		# If we filled missing fields, write back
		if updated:
			print(f"Updating {path_record} with new skymap info...")
			os.makedirs(path_output, exist_ok=True)
			with open(path_record, 'w') as f:
				json.dump(record, f, indent=2)
		else:
			print(f"No updates needed for {path_record}")

		# %%
		#	Initialize the variables for trigger
		gecko_digestor_trigger = False
		em_counterpart = False
		confidence_limit = 0.9 # 0-1 [%]
		obs_request_day = 0
		# gecko_priority: int = 0

		# %% [markdown]
		# - Alert type
		# 	- EARLYWARNING,PRELIMINARY,INITIAL,UPDATE,RETRACTION
		# - Significant
		# 	- true if trials factor Ã FAR < 1/month for CBC events, otherwise false
		# 	- true if trials factor Ã FAR < 1/year for burst events, otherwise false

		# %%
		# Determine event type by superevent_id prefix
		superevent_id = record['superevent_id']
		alert_type = record['alert_type']

		if superevent_id.startswith("S"):
			print("This is a superevent")
			slack = True
			superevent = True
			gecko_digestor_trigger = True

		elif superevent_id.startswith("MS"):
			# MS-prefix events are not superevents
			print("This is not a superevent")
			slack = False
			superevent = False
			gecko_digestor_trigger = False

		else:
			# Unknown prefix
			print("I don't know what this is")
			slack = True
			superevent = False
			gecko_digestor_trigger = False
		
		print(f"EVENT ID  : {superevent_id}")
		print(f"ALERT TYPE: {alert_type}")
		print(f"SLACK     : {slack}")
		print(f"PROCESS   : {gecko_digestor_trigger}")
		
		#	Alert Type Classification
		if (gecko_digestor_trigger) & (alert_type in ['EARLYWARNING','PRELIMINARY','INITIAL','UPDATE','RETRACTION']):
			#============================================================
			#	Designated alert type
			#============================================================
			#	No GECKO Digestor Trigger
			#------------------------------------------------------------
			#	Early warning
			if alert_type == 'EARLYWARNING':
				gecko_digestor_trigger = False
				# confidence_limit = 0.9 # 0-1 [%]
				# obs_request_day = 1
				print("This is an EARLYWARNING alert, stay tuned!")
			#	Preliminary
			elif alert_type == 'PRELIMINARY':
				gecko_digestor_trigger = True
				confidence_limit = 0.9 # 0-1 [%]
				# obs_request_day = 1
			#    Initial
			elif alert_type == 'INITIAL':
				gecko_digestor_trigger = True
				confidence_limit = 0.9 # 0-1 [%]
				# obs_request_day = 1
			#    Update
			elif alert_type == 'UPDATE':
				gecko_digestor_trigger = True
				confidence_limit = 0.9 # 0-1 [%]
				# obs_request_day = 3
			#	Retraction
			elif alert_type == 'RETRACTION':
				print("This is a RETRACTION alert, OH COME ON!")
		# else:
		# 	print("This is a NON-SIGNIFICANT alert")
		print(f"--> GECKO Digestor Trigger: {gecko_digestor_trigger}")
		print(f"--> Confidence limit      : {confidence_limit:.1%}")
		# %%
		if gecko_digestor_trigger:
			#============================================================
			#	Light Curve
			#============================================================
			#	Phase Calculation
			#------------------------------------------------------------
			print(f"Calculate SUN position...")
			t_gw = Time(record['event']['time'])
			t_now = Time.now()

			# 주어진 UTC 시간
			utc_time = Time(t_now.utc, scale='utc')

			# 태양의 적경, 적위 계산
			sun_position = get_sun(utc_time)
			sun_ra = sun_position.ra
			sun_dec = sun_position.dec

			# 결과 출력
			#
			if debug_mode == True:
				print("This is a DEBUG mode")
				phase = 1.0
			else:
				phase = t_now.jd - t_gw.jd
			print(f"Now     : {t_now.jd:.1f} JD")
			print(f"GW event: {t_gw.jd:.1f} JD")
			print("Sun: (RA, Dec):", sun_ra, sun_dec)
			print(f"Phase: {phase:.1f} days")
			# %%
			if record['event']['group'] == 'CBC':
				most_probable_event_prob = record['event']['classification'][most_probable_event]
				#   Yes EM counterpart
				if most_probable_event in ['BNS', 'NSBH',]:
					#   Possibilities to have an emission
					if (record['event']['properties']['HasNS']>0.9) & (record['event']['properties']['HasRemnant']>0.9):
						em_counterpart = True
				elif most_probable_event in ['BBH']:
					pass
				elif most_probable_event in ['Terrestrial']:
					gecko_digestor_trigger = False
					slack = False
					pass
			elif record['event']['group'] == 'Burst':
				most_probable_event = 'Burst'
				most_probable_event_prob = 1.0
				gecko_digestor_trigger = True
				slack = False 
			
			print(f"Most probable event: {most_probable_event} ({most_probable_event_prob*1e2:.1f}%)")

			# %%
			if record['external_coinc'] != None:
				if 'combined_skymap' in list(record['external_coinc'].keys()):
					skymap_str = record['external_coinc']['combined_skymap']
					skymap = read_skymap_bytes_to_table(skymap_str)
					print(f"{record['external_coinc']['search']} by {record['external_coinc']['observatory']} (t=t0+{record['external_coinc']['time_difference']} sec)")
			else:
				print(f"No external_coinc")

			# %%
			print(f"Initializing skymap RA and DEC...")
			skymap['RA'] = 0.0
			skymap['DEC'] = 0.0
			skymap['nside'] = np.zeros(len(skymap), dtype=np.int64)

			for i in np.arange(len(skymap)):
				uniq = skymap[i]['UNIQ']
				level, ipix = ah.uniq_to_level_ipix(uniq)
				nside = ah.level_to_nside(level)
				ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')

				skymap['RA'][i] = ra.deg
				skymap['DEC'][i] = dec.deg
				skymap['nside'][i] = nside

			# - Area within confidence region
			print(f"Calculate area within confidence region (50%, 90%)...")
			skymap.sort('PROBDENSITY', reverse=True)
			level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
			pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
			prob = pixel_area * skymap['PROBDENSITY']
			cumprob = np.cumsum(prob)
			i = cumprob.searchsorted(confidence_limit)
			j = cumprob.searchsorted(0.5) # 50%
			area_90 = pixel_area[:i].sum()
			area_50 = pixel_area[:j].sum()
			#	New column about cumulative probability
			skymap['CUMPROBDENSITY'] = cumprob

			print(f"{area_90.to_value(u.deg**2):1.3f} deg2 ({confidence_limit*1e2:.1f} %)")
			print(f"{area_50.to_value(u.deg**2):1.3f} deg2 ({0.5*1e2:.1f} %)")

			# %%
			print(f"Calculate RA and DEC ranges...")
			ramin, ramax = np.min(skymap['RA'][cumprob<confidence_limit]), np.max(skymap['RA'][cumprob<confidence_limit])
			decmin, decmax = np.min(skymap['DEC'][cumprob<confidence_limit]), np.max(skymap['DEC'][cumprob<confidence_limit])

			print(f"RA ranges: {ramin:.3f}~{ramax:.3f} deg")
			print(f"DEC ranges: {decmin:.3f}~{decmax:.3f} deg")

			ra_hms, dec_dms = deg2hmsdms(ramax, decmax)
			# %% [markdown]
			# - Check distance and most probable sky location
			# - Meta data from skymap

			print(f"Calculate most probable sky location...")
			level, ipix = ah.uniq_to_level_ipix(
				skymap[np.argmax(skymap['PROBDENSITY'])]['UNIQ']
			)
			ra, dec = ah.healpix_to_lonlat(ipix, ah.level_to_nside(level),
											order='nested')
			print(f'Most probable sky location (RA, Dec) = ({ra.deg:.3f}, {dec.deg:.3f})')

			# Print some information from FITS header
			try:
				distmean, diststd = skymap.meta['DISTMEAN'], skymap.meta['DISTSTD']
			except KeyError as e:
				print(f"Error: {e}. Using default distance values from configuration (d={TEMP_DIST_MEAN:.1f} +/- {TEMP_DIST_STD:.1f})")
				distmean, diststd = TEMP_DIST_MEAN, TEMP_DIST_STD
			print(f'Distance = {distmean:.3f} +/- {diststd:.3f} Mpc')

			# %%
			print(f"Calculate expected apparent magnitude...")
			plt.close()
			fig = plt.figure(figsize=(8, 5))
			filterlist = ['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K',]
			expected_magdict = expect_AT2017gfo(dprime=distmean, phase=phase, path_template_lc=TEMPLATE_LC, filterlist=filterlist, plot=True)
			yu, yl = plt.ylim()

			#	SMNet
			# for exptime, ls in zip([10, 60], ['--', '-']):
			# 	depth = np.median(gcktbl[f'depth_{exptime}min'][indx_smnet])
			# 	if (depth > yl) & (depth < yu):
			# 		plt.axhline(y=depth, ls=ls, lw=2, color='k', zorder=0)
			# 		plt.text(0.1, depth-0.1, f'1-m Tel. ({exptime}m)')
			# plt.axhline(y=20, ls=ls, lw=2, color='k', zorder=999)
			# plt.text(0.1, 20-0.1, f'1-m Tel.')

			#	KMTNet
			# for exptime, ls in zip([4, 24], ['--', '-']):
			# 	depth = np.median(gcktbl[f'depth_{exptime}min'][gcktbl['obs']=='KMTNet'])
			# 	# print(exptime, depth)
			# 	if (depth > yl) & (depth < yu):
			# 		plt.axhline(y=depth, ls=ls, lw=2, color='dodgerblue', zorder=0)
			# 		plt.text(0.1, depth-0.1, f'KMTNet ({exptime}m)')
			# plt.axhline(y=23, ls=ls, lw=2, color='dodgerblue', zorder=999)
			# plt.text(0.1, 23-0.1, f'KMTNet (10min)')

			lcpng = f"{path_output}/lc.png"
			print(f"Save KN light curve: {lcpng}")
			plt.savefig(lcpng)
			plt.close()

			# %%
			# - Read GLADE+ catalog (.fits for the faster I/O)

			'''- 1	GLADE no	GLADE+ catalog number
			- 2	PGC no	Principal Galaxies Catalogue number
			- 3	GWGC name	Name in the GWGC catalog
			- 4	HyperLEDA name	Name in the HyperLEDA catalog
			- 5	2MASS name	Name in the 2MASS XSC catalog
			- 6	WISExSCOS name	Name in the WISExSuperCOSMOS catalog (wiseX)
			- 7	SDSS-DR16Q name	Name in the SDSS-DR16Q catalog
			- 8	Object type flag
				- Q: the source is from the SDSS-DR16Q catalog
				- G:the source is from another catalog and has not been identified as a quasar
			- 9	RA	Right ascension in degrees
			- 10	Dec	Declination in degrees
			- 11	B	Apparent B magnitude
			- 12	B_err	Absolute error of apparent B magnitude
			- 13	B flag
				- 0: the B magnitude is measured
				- 1: the B magnitude is calculated from the B_J magnitude
			- 14	B_Abs	Absolute B magnitude
			- 15	J	Apparent J magnitude
			- 16	J_err	Absolute error of apparent J magnitude
			- 17	H	Apparent H magnitude
			- 18	H_err	Absolute error of apparent H magnitude
			- 19	K	Apparent K_s magnitude
			- 20	K_err	Absolute error of apparent K_s magnitude
			- 21	W1	Apparent W1 magnitude
			- 22	W1_err	Absolute error of apparent W1 magnitude
			- 23	W2	Apparent W2 magnitude
			- 24	W2_err	Absolute error of apparent W2 magnitude
			- 25	W1 flag
				- 0: the W1 magnitude is measured
				- 1: the W1 magnitude is calculated from the K_s magnitude
			- 26	B_J	Apparent B_J magnitude
			- 27	B_J err	Absolute error of apparent B_J magnitude
			- 28	z_helio	Redshift in the heliocentric frame
			- 29	z_cmb	Redshift converted to the Cosmic Microwave Background (CMB) frame
			- 30	z flag
				- 0: the CMB frame redshift and luminosity distance values given in columns 29 and 33 are not corrected for the peculiar velocity
				- 1: they are corrected values
			- 31	v_err	Error of redshift from the peculiar velocity estimation
			- 32	z_err	Measurement error of heliocentric redshift
			- 33	d_L	Luminosity distance in Mpc units
			- 34	d_L err	Error of luminosity distance in Mpc units
			- 35	dist flag
				- 0: the galaxy has no measured redshift or distance value
				- 1: it has a measured photometric redshift from which we have calculated its luminosity distance
				- 2: it has a measured luminosity distance value from which we have calculated its redshift
				- 3: it has a measured spectroscopic redshift from which we have calculated its luminosity distance
			- 36	M*	Stellar mass in 10^10 M_Sun units
			- 37	M*_err	Absolute error of stellar mass in 10^10 M_Sun units
			- 38	M* flag
				- 0: if the stellar mass was calculated assuming no active star formation
				- 1: if the stellar mass was calculated assuming active star formation
			- 39	Merger rate	Base-10 logarithm of estimated BNS merger rate in the galaxy in Gyr^-1 units
			- 40	Merger rate error	Absolute error of estimated BNS merger rate in the galaxy'''


			# %%
			print(f"Crossmatch Skymap and GLADE+ catalog...")
			result = crossmatch(skymap, coordinates, contours=(0.5, 0.9))
			if np.isnan(result.probdensity_vol).any() & np.isnan(result.searched_prob_vol).any():
				
				f = open(f"{path_output}/note.txt", "w")
				f.write(f"Crossmatch failed for P_3D. Use P_2D")
				f.close()

				print(f"Crossmatch failed for P_3D. Use P_2D")
				cat['searched_prob_vol'] = result.searched_prob
				cat['prob_vol'] = result.probdensity
			else:
				cat['searched_prob_vol'] = result.searched_prob_vol
				cat['prob_vol'] = result.probdensity_vol

			# %%
			print(f"Select galaxies within 90% confidence limit...")
			indx_vol90 = np.where(cat['searched_prob_vol']<confidence_limit)
			select_cat = cat[indx_vol90]
			#	Check bad matching
			print(f"Check bad matching...")
			if len(select_cat) == len(cat):
				print(f"Detected bad matching.")
				print(f"{len(select_cat)} galaxies matched (whole GLADE+ catalog)")
				indx_tmp = np.where(
					# (select_cat['col33']<distmean+diststd) &
					# (select_cat['col33']>distmean-diststd) &
					(select_cat['col9']<ramax) &
					(select_cat['col9']>ramin) &
					(select_cat['col10']<decmax) &
					(select_cat['col10']>decmin)
					)
				select_cat = select_cat[indx_tmp]
				if os.path.exists(f"{path_output}/note.txt"):
					#	add
					write_mode = "a"
				else:
					#	write new one
					write_mode = "w"
				f = open(f"{path_output}/note.txt", write_mode)
				f.write(f"\nDetected bad matching\n")
				f.close()

			#	Stellar mass
			select_cat['stellar_mass'] = select_cat['col36']
			select_cat['flag_stmass'] = ~select_cat['col36'].mask

			#	Check size of masked rows
			select_cat['stellar_mass'][select_cat['col36'].mask] = np.min(select_cat['col36'])

			#	Extra probability
			_prob_value = select_cat['stellar_mass']*select_cat['prob_vol']
			select_cat['prob_vol_x_stmass'] = (_prob_value-np.min(_prob_value))/(np.max(_prob_value)-np.min(_prob_value))
			print(f"{len(select_cat):,} galaxies selected")

			# %%
			# - LIGHT version of Catalog
			#	Key to sort the candidate
			probkey = "prob_vol"
			# probkey = "prob_vol_x_stmass"

			print(f"="*60)
			print(f"Generate Light Version of Host Galaxy Candidate Catalog")
			print(f"="*60)
			simple_galcat = Table()
			simple_galcat['name'] = select_cat['col1']
			simple_galcat['ra'] = select_cat['col9']
			simple_galcat['dec'] = select_cat['col10']
			simple_galcat['d_L'] = select_cat['col33']
			simple_galcat['prob_vol'] = select_cat['prob_vol']
			simple_galcat['stmass'] = select_cat['stellar_mass']
			simple_galcat['prob_vol_x_stmass'] = select_cat['prob_vol_x_stmass']

			#	To prevent all values from being NaN.
			print(f"Preventing all values from being NaN...")
			simple_galcat['prob_vol_x_stmass'] = np.nan_to_num(
				simple_galcat['prob_vol_x_stmass'],
				nan=np.min(simple_galcat['prob_vol_x_stmass'][(simple_galcat['prob_vol_x_stmass']!=0.0) & ~np.isnan(simple_galcat['prob_vol_x_stmass'])])
				)

			#	Sort by prob --> rank
			print(f"Sorting by {probkey}...")
			simple_galcat = simple_galcat[np.flipud(np.argsort(simple_galcat[probkey]))]
			simple_galcat['rank'] = np.arange(len(simple_galcat), dtype=int)

			#   Formatting
			print(f"Formatting...")
			for key in simple_galcat.keys():
				# if key in ['ra', 'dec', 'd_L', probkey]:
				if key in ['ra', 'dec', 'd_L',]:
					simple_galcat[key].format = '.3f'

			#   Cumulative Probability
			print(f"Calculating cumulative probability...")
			cumsum_prob_gal = np.cumsum(simple_galcat[probkey])
			sum_prob_gal = np.sum(simple_galcat[probkey])

			simple_galcat['confidence'] = 0.0
			simple_galcat['confidence'][np.max(cumsum_prob_gal)*1.0>=cumsum_prob_gal] = 1.0
			simple_galcat['confidence'][np.max(cumsum_prob_gal)*0.99>=cumsum_prob_gal] = 0.99
			simple_galcat['confidence'][np.max(cumsum_prob_gal)*0.95>=cumsum_prob_gal] = 0.95
			simple_galcat['confidence'][np.max(cumsum_prob_gal)*0.9>=cumsum_prob_gal] = 0.9
			simple_galcat['confidence'][np.max(cumsum_prob_gal)*0.5>=cumsum_prob_gal] = 0.5

			simple_galcat['obj'] = record['superevent_id']
			simple_galcat['note'] = simple_galcat['name']
			simple_galcat['weight'] = simple_galcat[probkey]

			#   Meta data
			print(f"Adding meta data...")
			galcat_meta_dict = {
				'superevent_id': record['superevent_id'],
				'alert_type': record['alert_type'],
				'most_probable_event': most_probable_event,
				'most_probable_event_prob': most_probable_event_prob,
				'confidence': confidence_limit,
				'ordering': probkey
			}
			simple_galcat.meta = galcat_meta_dict

			print(f"Saving the host galaxy candidate catalog...")
			simple_galcat_name = f"{path_output}/HostGalaxyCatalog_{confidence_limit*1e2:g}.csv"
			print(f"{simple_galcat_name}")
			simple_galcat.write(simple_galcat_name, format='csv', overwrite=True)

			#	ds9 region file
			simple_galcat_region = simple_galcat_name.replace("csv", "reg")
			print(f"Creating ds9 region file...")
			print(f"{simple_galcat_region}")
			create_ds9_regions(simple_galcat["ra"], simple_galcat["dec"], simple_galcat["name"], filename=simple_galcat_region)

			# %%
			# Prepare the healpix map (unchanged)
			# nside   = 2**7 # This fixed valude is only for the visualization
			# nside   = np.max(skymap['nside'])
			nside   = int(64)
			# nside   = skymap['nside']
			npix    = hp.nside2npix(nside)
			hp_map  = np.zeros(npix)
			ra_rad  = np.radians(skymap['RA'])
			dec_rad = np.radians(skymap['DEC'])
			indices = hp.ang2pix(nside, theta=0.5*np.pi - dec_rad, phi=ra_rad)
			hp_map[indices] = skymap['PROBDENSITY']

			cmap = cmap

			# Tick positions in radians for every 30°
			lon_ticks  = np.radians(np.arange(-150, 181, 30))
			lat_ticks  = np.radians(np.arange(-90,   91,  30))
			lon_labels = [f"{int(l)}°" for l in np.arange(-150, 181, 30)]
			lat_labels = [f"{int(l)}°" for l in np.arange(-90,   91,  30)]

			hp_mollview_title = f"{record['superevent_id']}_{record['alert_type']}: {most_probable_event} ({most_probable_event_prob:.1%})\n50% area: {area_50.to(u.deg**2):.1f}\n90% area: {area_90.to(u.deg**2):.1f}\nd={distmean:.0f}+{diststd:.0f} Mpc"

			# --- Figure 1: no rotation ---
			print(f"Drawing no rotation skymap...")
			fig1 = plt.figure(figsize=(10, 6))

			# hp.mollview(
			# 	hp_map,
			# 	# fig=fig1.number, sub=(1,1,1),
			# 	title=f"{record['superevent_id']}_{record['alert_type']} (no rot)",
			# 	cmap=cmap, cbar=True, flip='astro', unit=r'$\rm P_{3D}$'
			# )
			projview(
				hp_map, 
				coord=["E"], 
				#
				graticule=True, 
				graticule_labels=True, 
				#
				latitude_grid_spacing=30,
				longitude_grid_spacing=45,
				#
				cmap=cmap, 
				# cbar=True, 
				cbar=False, 
				title=hp_mollview_title,
				projection_type="mollweide",
				#
				custom_xtick_labels=["21h", "18h", "15h", "12h", "9h", "6h", "3h",],
				# custom_ytick_labels=[""],
				#
				flip='astro', 
				unit=r'$\rm P_{3D}$',
				rot=(180, 0),
				# rot_graticule=True,

			)
			# hp.graticule(dpar=30, dmer=30, coord='E', c='k', alpha=0.75)

			# 2) 사각형 네 꼭짓점 (deg 단위)
			lon1, lon2 =  1,  2  # 예시 경도 (°)
			lat1, lat2 = -1,  1  # 예시 위도 (°)

			# 마지막에 첫 점을 다시 넣어 닫기
			lons = [lon1, lon1, lon2, lon2, lon1]
			lats = [lat1, lat2, lat2, lat1, lat1]

			# ax = plt.gca()
			# ax.set_xticks(lon_ticks)
			# ax.set_xticklabels(lon_labels)
			# ax.set_yticks(lat_ticks)
			# ax.set_yticklabels(lat_labels)

			plt.tight_layout()
			plt.savefig(f"{path_output}/skymap_no_rot.png", dpi=100)
			plt.close(fig1)
			# %%

			# --- Figure 2: rotated 90° ---
			print(f"Drawing rotated skymap...")
			fig2 = plt.figure(figsize=(10, 4))

			hp.mollview(
				hp_map,
				# fig=fig2.number, sub=(1,1,1),
				# title=f"{record['superevent_id']}_{record['alert_type']} (rot 90°)",
				cmap=cmap, cbar=True, flip='astro', unit=r'$\rm P_{3D}$',
				title=hp_mollview_title,
				rot=(180, -45), 
			)
			# projview(
			# 	hp_map, coord=["E"], graticule=True, graticule_labels=True, projection_type="mollweide",
			# 	title=f"{record['superevent_id']}_{record['alert_type']}: 50% (90%) area: {area_50.to(u.deg**2):.1f} ({area_90.to(u.deg**2):.1f})",
			# 	cmap=cmap, cbar=True, flip='astro', unit=r'$\rm P_{3D}$',
			# 	rot=(180, -45)
			# )
			hp.graticule(dpar=30, dmer=30, coord='E', c='k', alpha=0.75)

			ax = plt.gca()


			plt.tight_layout()
			plt.savefig(f"{path_output}/skymap_rot90.png", dpi=100)
			plt.close(fig2)			

			# %%
			plt.close()
			print(f"Drawing cumulative probability...")
			fig = plt.figure(figsize=(8, 6))
			norm_factor = np.max(cumsum_prob_gal)
			plt.plot(cumsum_prob_gal/norm_factor, '-', mfc='w', mew=3, ms=10, lw=3, c='g', label=f"All ({len(cumsum_prob_gal)})", zorder=0)
			plt.axhline(y=0.99*sum_prob_gal/norm_factor, ls='-', c='tomato', lw=3, alpha=1.0, label=f"99% ({len(simple_galcat[simple_galcat['confidence']<=0.99]):,})")
			plt.axhline(y=0.95*sum_prob_gal/norm_factor, ls='-.', c='tomato', lw=3, alpha=0.75, label=f"95% ({len(simple_galcat[simple_galcat['confidence']<=0.95]):,})")
			plt.axhline(y=0.9*sum_prob_gal/norm_factor, ls=':', c='tomato', lw=3, alpha=0.75, label=f"90% ({len(simple_galcat[simple_galcat['confidence']<=0.9]):,})")
			plt.axhline(y=0.5*sum_prob_gal/norm_factor, ls='--', c='tomato', lw=3, alpha=0.5, label=f"50% ({len(simple_galcat[simple_galcat['confidence']<=0.5]):,})")
			plt.xlabel('Number of Host Candidates', fontsize=20)

			if probkey == 'prob_vol':
				ylabel = r'Cumulative $\rm P_{3D}$'
			elif probkey == 'prob_vol_x_stmass':
				ylabel = r'Cumulative $\rm P_{3D}xP_{M_{*}}$'
			plt.ylabel(ylabel)
			# plt.ylabel(r'Cumulative $\rm P_{3D}$', fontsize=20)
			yl, yu = plt.ylim()
			plt.ylim([0, yu])
			plt.legend(loc='lower right', fontsize=20)
			# _ = plt.yticks(np.arange(0, 1.1, 0.1))
			plt.grid('both', color='silver', ls='--', lw=1, alpha=0.5)
			yticks = plt.yticks(fontsize=20)
			plt.xticks(fontsize=20)
			plt.yticks(np.arange(0, 1.01, 0.1), fontsize=20)

			plt.tight_layout()
			plt.savefig(f"{path_output}/cumulative_p3d_HostGalaxy.png", dpi=100,)
			plt.close(fig)

			# %%
			print(f"-"*60)
			print("Setting Tiling Pattern")
			print(f"-"*60)
			#	Default values for escaping the error

			probkey = "sum_p3d"
			n_top = 300 # number of tiles
			nside   = ah.level_to_nside(level)

			tile_number_dict = {}
			for obs in telescopes4tile:
				print(f"Generating Tiling Pattern for {obs}")
				# 1) load & filter
				print(f"Loading skygrid...")
				skygrid_cat = load_skygrid(obs)
				skygrid_cat = clean_skygrid(skygrid_cat)
				print(f"Filtering skygrid...")
				# selected_tiles = process_tiles(skygrid_cat, skymap, simple_galcat, confidence_limit, nside)

				idx, sums = select_tiles_and_sum_prob(
					skygrid_cat,
					skymap,
					confidence_limit=confidence_limit,
					nside_high=np.max(nside),
				)

				# 선택된 타일 테이블만 떼어내고
				selected_tiles = skygrid_cat[idx]
				# 그리고 컬럼으로 합산 PROBDENSITY를 붙여주면 끝
				selected_tiles['sum_p3d'] = sums
				selected_tiles['nor_sum_p3d'] = sums/np.max(sums)

				selected_tiles = selected_tiles[np.argsort(-1*selected_tiles['nor_sum_p3d'])]

				# if len(selected_tiles)==0:
					# continue
				# 3) save & plot
				# 1) CSV 저장
				print(f"Saving selected skygrid catalog...")
				save_skygrid_catalog(selected_tiles, obs, record, probkey, confidence_limit, path_output)

				# 2) Tiling map 그림
				print(f"Drawing tiling map...")
				plot_tiling_map(selected_tiles[:n_top], hp_map, obs, record, path_output)

				# 2) Tiling map zoom 그림
				print(f"Drawing tiling map zoom...")
				plot_tiling_map_zoom(
					sel=selected_tiles[:n_top],
					skymap=hp_map,
					obs=obs,
					record=record,
					path_output=path_output,
					ra_center=record['ramax'],
					dec_center=record['decmax'],
					fov_deg=3000,
					xsize=300,
				)
				# 3) Cumulative probability 그림
				print(f"Drawing cumulative probability...")
				plot_cumulative_dist(selected_tiles, obs, probkey, path_output)

				print(f"Number of selected tiles: {len(selected_tiles)} ({confidence_limit:.0%})")
				tile_number_dict[obs] = len(selected_tiles)
				# break
			# %%
			#============================================================
			#   Summary
			#------------------------------------------------------------
			delt = time.time() - st

			import textwrap

			summary_txt = textwrap.dedent(f"""\
			Event Information:
			Event ID: {record['superevent_id']}
			Alert Type: {record['alert_type']}
			Trigger Time: {record['event']['time']}
			Processing Start: {t_now.isot}Z
			Processing Time: {delt:.1f} sec

			Localization:
			Phase: {phase:1.1f} days
			Classification: {most_probable_event} ({most_probable_event_prob*100:.1f}%)
			Distance: {distmean:.1f} +/- {diststd:.1f} Mpc

			Sky Position:
			RA, Dec (deg): {ramax:1.3f},{decmax:1.3f}
			RA, Dec (hms, dms): {ra_hms},{dec_dms}

			Sky Area:
			90% Area: {area_90.to(u.deg**2).value:.1f} deg²
			50% Area: {area_50.to(u.deg**2).value:.1f} deg²

			Host Galaxies:
			Total: {len(simple_galcat):,}
			Within 90%: {len(simple_galcat[simple_galcat['confidence']<=0.9]):,}
			Within 50%: {len(simple_galcat[simple_galcat['confidence']<=0.5]):,}

			Expected Magnitudes:
			""")

			# Add expected magnitudes
			for filte in list(expected_magdict.keys()):
				summary_txt += f"{filte} band: {expected_magdict[filte]:.1f} mag\n"

			summary_txt += "\nObservation Coverage:\n"
			for telescope, count in tile_number_dict.items():
				summary_txt += f"{telescope}: {count:,} tiles\n"
			
			summary_txt += f"\nProcessed Time: {delt/60.:.1f} min\n"
			summary_txt += f"\nOutput Path: {path_output}\n"

			# Format and write summary
			summary_file = f"{path_output}/summary.txt"
			with open(summary_file, "w") as f:
				f.write(summary_txt)

			print(summary_txt)
			# if os.path.exists(summary_file):
			# 	os.system(f"cat {summary_file}")
			# else:
			# 	print(f"The Summary File {summary_file} does not exist")
			
			# # Fin. ALERT!!!
			# %%
			# 0. debug 모드에서는 무조건 슬랙 전송 비활성화
			if debug_mode:
				slack = False

			# 1. slack=True 이면서 debug_mode=False 인 경우만 여기로 들어옵니다.
			if slack:
				# (a) 이미지 파일 경로 설정
				img_file = f"{path_output}/skymap_no_rot.png"

				# (b) Process done 메시지+이미지 업로드
				initial_comment = (
					f"[`GeckoDigestor`] {record['superevent_id']}-{record['alert_type']}, "
					f"({most_probable_event} {most_probable_event_prob:.1%}"
				)
				if most_probable_event != 'Burst':
					initial_comment += f", d={distmean:.1f}±{diststd:.1f} Mpc"
				initial_comment += ")"

				resp = client.files_upload(
					channels=slack_channel,
					file=img_file,
					initial_comment=initial_comment
				)

				# (c) 업로드된 메시지에서 thread_ts 추출
				shares     = resp['file'].get('shares', {}).get('public', {})
				channel_id = next(iter(shares.keys()))
				thread_ts  = shares[channel_id][0]['ts']

				# (d) 같은 스레드에 summary 텍스트 포스팅
				client.chat_postMessage(
					channel=channel_id,
					thread_ts=thread_ts,
					text=f"```{summary_txt}```"
				)

				print(f"DONE! ({delt/60.:.1f} min)")

			# 2. slack=False 이거나 debug_mode=True → 전혀 보내지 않음
			else:
				print(f"[slack={slack}] 슬랙 메시지를 보내지 않습니다.")
				print("Waiting for a GW alert...")
				
		#	Update the event log
		print(f"Updating the event log...")
		indx_match = eventlogtbl['output_dir'] == os.path.basename(path_output)
		eventlogtbl['processed'][indx_match] = True

		# eventlogtbl.add_row(
		# 	[record['superevent_id'], record['alert_type'], most_probable_event, record['ramax'], record['decmax'], record['area_90'], record['distmean'], record['diststd'], path_output, True]
		# )
		# for key in eventlogtbl.keys():
		# 	if key in ['ramax', 'decmax', 'area_90', 'distmean', 'diststd']:
		# 		eventlogtbl[key].format = ".3f"

		# print(eventlogtbl[-1])
		eventlogtbl.write(f"{path_out}/event.log", format='ascii.fixed_width', overwrite=True)

		print(f"DONE!")

	time.sleep(TIME_TO_SLEEP)
