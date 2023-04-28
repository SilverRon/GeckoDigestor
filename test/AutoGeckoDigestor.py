# %% [markdown]
# # Library

# %%
import time

import astropy_healpix as ah
import numpy as np
from astropy.table import Table

from astropy.coordinates import SkyCoord
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap.postprocess import crossmatch
from astropy import units as u
from astropy.io import ascii
from astropy.time import Time

# %%
#	Append the path to the GECKO project root directory
import sys
sys.path.append('..')
from config.config import *
from src.AlertReceiver import *
from src.GeckoHelper import *
from util.util import *
#%%
#	Scheduler
sys.path.append('../observation')
sys.path.append('../observation/Obsscheduler')
from Obsscheduler.obsscheduler import ObsScheduler


# %%
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

# %% [markdown]
# # 1. Receive the GW alert

# %% [markdown]
# - Online Test
path_out = '../output'
slack = True
#	Read HEAVY GLADE+ catalog in advance
try:
	cat
except NameError:
	cat = Table.read('../data/GLADE+_230419.fits', format='fits',)				
coordinates = SkyCoord(cat['col9']*u.deg, cat['col10']*u.deg, cat['col33']*u.Mpc)

#	GECKO facilities
gcktbl = ascii.read("../data/gecko.facilities.tsv")
gcktbl['depth_3min'] = scale_depth(gcktbl['depth'], gcktbl['exptime'], 3.)
gcktbl['depth_4min'] = scale_depth(gcktbl['depth'], gcktbl['exptime'], 4.)
gcktbl['depth_10min'] = scale_depth(gcktbl['depth'], gcktbl['exptime'], 10.)
gcktbl['depth_24min'] = scale_depth(gcktbl['depth'], gcktbl['exptime'], 24.)
gcktbl['depth_30min'] = scale_depth(gcktbl['depth'], gcktbl['exptime'], 30.)
gcktbl['depth_60min'] = scale_depth(gcktbl['depth'], gcktbl['exptime'], 60.)
indx_smnet = np.where(
    (gcktbl['obs']=='LOAO') |
    (gcktbl['obs']=='CBNUO') |
    (gcktbl['obs']=='SAO') |
    (gcktbl['obs']=='KCT') |
    (gcktbl['obs']=='RASA36') |
    (gcktbl['obs']=='DOAO') |
    (gcktbl['obs']=='KHAO') |
    (gcktbl['obs']=='MDFTS') |
    (gcktbl['obs']=='LSGT') |
    (gcktbl['obs']=='MAAO') |
    (gcktbl['obs']=='DNSM')
)
# %%
#	Online Test
config = {
	'group.id': '',
	# 'auto.offset.reset': 'earliest'
	}

consumer = Consumer(
	config=config,
	client_id=GCN_KAFKA_CONFIG['client_id'],
	client_secret=GCN_KAFKA_CONFIG['client_secret'],
	)
consumer.subscribe(['igwn.gwalert'])

while True:
	print("Waiting for a GW alert...")
	for message in consumer.consume(timeout=30*24*60*60):
		record, skymap = AlertReceiver(message.value())
		st = time.time()

		# 새로 생성할 디렉토리 이름
		new_directory = f"{path_out}/{record['superevent_id']}_{record['alert_type']}"
		# 디렉토리가 존재하지 않으면 생성
		if not os.path.exists(new_directory):
			os.makedirs(new_directory)
		# 딕셔너리를 JSON 문자열로 변환
		record_str = json.dumps(record)
		# JSON 문자열을 파일로 저장
		with open(f'{new_directory}/record.json', 'w') as file:
			file.write(record_str)

		# if (record['alert_type'] == 'RETRACTION') & (skymap == None):
		if (record['alert_type'] != 'RETRACTION'):
			write_skymap_to_fits(skymap, path_output=f"{new_directory}/skymap.fits")


			# %% [markdown]
			# - Manual Test

			# %%
			#	Offline Test
			# with open('../data/MS181101ab-preliminary.json', 'r') as f:
			# with open('../data/MS181101ab-earlywarning.json', 'r') as f:
			# with open('../data/MS181101ab-ext-update.json', 'r') as f:
			# with open('../data/MS181101ab-initial.json', 'r') as f:
			#     record = f.read()

			# record, skymap = AlertReceiver(record)

			# %% [markdown]
			# - Path

			# %%
			path_output = f"../output/{record['superevent_id']}_{record['alert_type']}"

			# %% [markdown]
			# # 2. Check a Singnificance

			# %% [markdown]
			# - Initialize the variables for trigger

			# %%
			gecko_digestor_trigger = False
			em_counterpart = False
			confidence_limit = 0.0 # 0-1 [%]
			obs_request_day = 0
			# gecko_priority: int = 0

			# %% [markdown]
			# - Alert type
			# 	- EARLYWARNING,PRELIMINARY,INITIAL,UPDATE,RETRACTION
			# - Significant
			# 	- true if trials factor Ã FAR < 1/month for CBC events, otherwise false
			# 	- true if trials factor Ã FAR < 1/year for burst events, otherwise false

			# %%
			# if record['event']['significant']:
			if record['alert_type'] in ['EARLYWARNING','PRELIMINARY','INITIAL','UPDATE','RETRACTION']:
				#============================================================
				#	Designated alert type
				#============================================================
				#	No GECKO Digestor Trigger
				#------------------------------------------------------------
				#	Early warning
				if record['alert_type'] == 'EARLYWARNING':
					gecko_digestor_trigger = False
					# confidence_limit = 0.9 # 0-1 [%]
					obs_request_day = 1
					print("This is an EARLYWARNING alert, stay tuned!")
				#	Preliminary
				elif record['alert_type'] == 'PRELIMINARY':
					gecko_digestor_trigger = True
					confidence_limit = 0.9 # 0-1 [%]
					obs_request_day = 1
				#    Initial
				elif record['alert_type'] == 'INITIAL':
					gecko_digestor_trigger = True
					confidence_limit = 0.9 # 0-1 [%]
					obs_request_day = 1
				#    Update
				elif record['alert_type'] == 'UPDATE':
					gecko_digestor_trigger = True
					confidence_limit = 0.9 # 0-1 [%]
					obs_request_day = 3
				#	Retraction
				elif record['alert_type'] == 'RETRACTION':
					print("This is a RETRACTION alert, OH COME ON!")
			# else:
			# 	print("This is a NON-SIGNIFICANT alert")

			print(f"{record['alert_type']} event detected")
			if gecko_digestor_trigger:
				print(f"--> Initiate GECKO Digestor Trigger")
				print(f"--> Set confidence limit to {confidence_limit*1e2:.1f}%")
				#============================================================
				#	Light Curve
				#============================================================
				#	Phase Calculation
				#------------------------------------------------------------
				t_gw = Time(record['event']['time'])
				# t_now = Time("2018-11-03T04:30:46.654Z")
				t_now = Time.now()
				phase = t_now.jd - t_gw.jd
				print(f"Phase: {phase:.1f} days")
				# %%
				most_probable_event = max(record['event']['classification'], key=record['event']['classification'].get)
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
					pass
				print(f"Most probable event: {most_probable_event} ({most_probable_event_prob*1e2:.1f}%)")
				print(f"Gecko Digestor Trigger: {gecko_digestor_trigger}")

				# %% [markdown]
				# - Check the existence of an `external_coinc`

				# %%
				if record['external_coinc'] != None:
					if 'combined_skymap' in list(record['external_coinc'].keys()):
						skymap_str = record['external_coinc']['combined_skymap']
						skymap = read_skymap_bytes_to_table(skymap_str)
						print(f"{record['external_coinc']['search']} by {record['external_coinc']['observatory']} (t=t0+{record['external_coinc']['time_difference']} sec)")
				else:
					print(f"No external_coinc")

				# %% [markdown]
				# ### Skymap Analysis

				# %%
				skymap['RA'] = 0.0
				skymap['DEC'] = 0.0

				for i in np.arange(len(skymap)):
					uniq = skymap[i]['UNIQ']
					level, ipix = ah.uniq_to_level_ipix(uniq)
					nside = ah.level_to_nside(level)
					ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')

					skymap['RA'][i] = ra.deg
					skymap['DEC'][i] = dec.deg

				# %% [markdown]
				# - Area within confidence region

				# %%
				skymap.sort('PROBDENSITY', reverse=True)
				level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
				pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
				prob = pixel_area * skymap['PROBDENSITY']
				cumprob = np.cumsum(prob)
				i = cumprob.searchsorted(confidence_limit)
				area_90 = pixel_area[:i].sum()
				#	New column about cumulative probability
				skymap['CUMPROBDENSITY'] = cumprob

				print(f"{area_90.to_value(u.deg**2):1.3f} deg2 ({confidence_limit*1e2:.1f} %)")

				# %%
				ramin, ramax = np.min(skymap['RA'][cumprob<confidence_limit]), np.max(skymap['RA'][cumprob<confidence_limit])
				decmin, decmax = np.min(skymap['DEC'][cumprob<confidence_limit]), np.max(skymap['DEC'][cumprob<confidence_limit])

				print(f"RA ranges: {ramin:.3f}~{ramax:.3f} deg")
				print(f"DEC ranges: {decmin:.3f}~{decmax:.3f} deg")

				# %% [markdown]
				# - Check distance and most probable sky location
				# - Meta data from skymap

				# %%
				level, ipix = ah.uniq_to_level_ipix(
					skymap[np.argmax(skymap['PROBDENSITY'])]['UNIQ']
				)
				ra, dec = ah.healpix_to_lonlat(ipix, ah.level_to_nside(level),
												order='nested')
				print(f'Most probable sky location (RA, Dec) = ({ra.deg:.3f}, {dec.deg:.3f})')

				# Print some information from FITS header
				distmean, diststd = skymap.meta['DISTMEAN'], skymap.meta['DISTSTD']
				print(f'Distance = {distmean:.3f} +/- {diststd:.3f} Mpc')

				plt.close()
				fig = plt.figure(figsize=(8, 6))
				filterlist = ['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K',]
				expected_magdict = expect_AT2017gfo(dprime=distmean, phase=phase, filterlist=filterlist, plot=True)
				yu, yl = plt.ylim()

				#	SMNet
				for exptime, ls in zip([10, 60], ['--', '-']):
					depth = np.median(gcktbl[f'depth_{exptime}min'][indx_smnet])
					if (depth > yl) & (depth < yu):
						plt.axhline(y=depth, ls=ls, lw=2, color='k', zorder=0)
						plt.text(0.1, depth-0.1, f'1-m Tel. ({exptime}m)')

				#	KMTNet
				for exptime, ls in zip([4, 24], ['--', '-']):
					depth = np.median(gcktbl[f'depth_{exptime}min'][gcktbl['obs']=='KMTNet'])
					# print(exptime, depth)
					if (depth > yl) & (depth < yu):
						plt.axhline(y=depth, ls=ls, lw=2, color='dodgerblue', zorder=0)
						plt.text(0.1, depth-0.1, f'KMTNet ({exptime}m)')

				plt.savefig(f"{path_output}/lc.png")
				plt.close()
				# %% [markdown]
				# # 3. Catalog Matching

				# %% [markdown]
				# - Read GLADE+ catalog (.fits for the faster I/O)

				# %% [markdown]
				# ### Read & Matching

				# %%
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
				indx_vol90 = np.where(cat['searched_prob_vol']<confidence_limit)
				select_cat = cat[indx_vol90]
				#	Stellar mass
				select_cat['stellar_mass'] = select_cat['col36']
				select_cat['flag_stmass'] = ~select_cat['col36'].mask
				select_cat['stellar_mass'][select_cat['col36'].mask] = np.min(select_cat['col36'])
				#	Extra probability
				select_cat['prob_vol_x_stmass'] = min_max_normalize(select_cat['stellar_mass']*select_cat['prob_vol'])

				# %%
				#	Key to sort the candidate
				# probkey = "prob_vol"
				probkey = "prob_vol_x_stmass"

				# %% [markdown]
				# - LIGHT version of Catalog

				# %%
				simple_galcat = Table()
				simple_galcat['name'] = select_cat['col1']
				simple_galcat['ra'] = select_cat['col9']
				simple_galcat['dec'] = select_cat['col10']
				simple_galcat['d_L'] = select_cat['col33']
				simple_galcat['prob_vol'] = select_cat['prob_vol']
				simple_galcat['stmass'] = select_cat['stellar_mass']
				simple_galcat[probkey] = select_cat[probkey]
				#	Sort by prob --> rank
				simple_galcat = simple_galcat[np.flipud(np.argsort(simple_galcat[probkey]))]
				simple_galcat['rank'] = np.arange(len(simple_galcat), dtype=int)
				#   Formatting
				for key in simple_galcat.keys():
					if key in ['ra', 'dec', 'd_L', probkey]:
						simple_galcat[key].format = '.3f'
				#   Cumulative Probability
				cumsum_prob_gal = np.cumsum(simple_galcat[probkey])
				sum_prob_gal = np.sum(simple_galcat[probkey])

				simple_galcat['confidence'] = 0.0
				simple_galcat['confidence'][np.max(cumsum_prob_gal)*1.0>=cumsum_prob_gal] = 1.0
				simple_galcat['confidence'][np.max(cumsum_prob_gal)*0.95>=cumsum_prob_gal] = 0.95
				simple_galcat['confidence'][np.max(cumsum_prob_gal)*0.9>=cumsum_prob_gal] = 0.9
				simple_galcat['confidence'][np.max(cumsum_prob_gal)*0.5>=cumsum_prob_gal] = 0.5

				simple_galcat['obj'] = record['superevent_id']
				simple_galcat['note'] = simple_galcat['name']
				simple_galcat['weight'] = simple_galcat[probkey]

				#   Meta data
				galcat_meta_dict = {
					'superevent_id': record['superevent_id'],
					'alert_type': record['alert_type'],
					'most_probable_event': most_probable_event,
					'most_probable_event_prob': most_probable_event_prob,
					'confidence': confidence_limit,
					'ordering': probkey
				}
				simple_galcat.meta = galcat_meta_dict

				simple_galcat_name = f"{path_output}/HostGalaxyCatalog_{confidence_limit*1e2:g}.csv"
				simple_galcat.write(simple_galcat_name, format='csv', overwrite=True)
				print(f"Save the host galaxy candidate catalog ({os.path.basename(simple_galcat_name)})")

				# %%
				#	Skymap
				# fig = plt.figure(figsize=(10, 8))
				fig = plt.figure(figsize=(20, 16))

				plt.scatter(skymap['RA'][cumprob<0.9], skymap['DEC'][cumprob<0.9], c=skymap['CUMPROBDENSITY'][cumprob<0.9], zorder=0, alpha=0.5)
				xl, xr = plt.xlim()
				cbar = plt.colorbar()
				cbar.set_label(r"$\rm P_{3D}$")

				plt.grid('both', ls='--', c='grey', alpha=0.5)
				plt.xlim([xr, xl])
				plt.xlabel('RA [deg]', fontsize=20)
				plt.ylabel('Dec [deg]', fontsize=20)
				plt.xticks(fontsize=20)
				plt.yticks(fontsize=20)
				plt.savefig(f"{path_output}/skymap_90.png", dpi=100,)
				plt.close()
				# %%
				# fig = plt.figure(figsize=(20, 8))
				fig = plt.figure(figsize=(40, 16))

				plt.subplot(121)

				plt.scatter(simple_galcat['ra'], simple_galcat['dec'], c=simple_galcat[probkey], cmap='hot', edgecolors='k')
				cbar_gal = plt.colorbar()

				plt.scatter(skymap['RA'][cumprob<0.9], skymap['DEC'][cumprob<0.9], c=skymap['CUMPROBDENSITY'][cumprob<0.9], zorder=0)
				xl, xr = plt.xlim()
				# cbar = plt.colorbar()

				plt.grid('both', ls='--', c='grey', alpha=0.5)
				plt.xlim([xr, xl])
				plt.xlabel('RA [deg]', fontsize=20)
				plt.ylabel('Dec [deg]', fontsize=20)
				plt.xticks(fontsize=20)
				plt.yticks(fontsize=20)
				# plt.tight_layout()

				plt.subplot(122)
				plt.plot(cumsum_prob_gal, '.-', mfc='w', mew=3, ms=10, lw=3, c='g', label=f"All ({len(cumsum_prob_gal)})", zorder=0)
				plt.axhline(y=0.95*sum_prob_gal, ls='-', c='tomato', lw=3, alpha=0.75, label=f"95% ({len(simple_galcat[simple_galcat['confidence']<=0.95])})")
				plt.axhline(y=0.9*sum_prob_gal, ls=':', c='tomato', lw=3, alpha=0.75, label=f"90% ({len(simple_galcat[simple_galcat['confidence']<=0.9])})")
				plt.axhline(y=0.5*sum_prob_gal, ls='--', c='tomato', lw=3, alpha=0.5, label=f"50% ({len(simple_galcat[simple_galcat['confidence']<=0.5])})")
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
				plt.xticks(fontsize=20)
				plt.yticks(fontsize=20)

				plt.tight_layout()
				plt.savefig(f"{path_output}/cumulative_p3d_HostGalaxy.png", dpi=100,)
				plt.close()
				#%%
				#============================================================
				#	Galaxy-targeted observation
				#------------------------------------------------------------
				#	LSGT, LOAO, KOREA FACILITIES, ...
				#============================================================
				date = Time.now()
				subsequent_days = obs_request_day
				file_prefix ='GECKO_'
				# for i in range(subsequent_days):
				for i in range(subsequent_days):
					#------------------------------------------------------------
					#	LSGT
					#------------------------------------------------------------
					obssch = ObsScheduler(
						target_db=simple_galcat[:NUMBER_GALAXY_LIMIT],
						name_telescope='LSGT',
						entire_night=True,
						duplicate_when_empty=False,
						date=date,
						autofocus_at_init=False
						)
					# obssch.write_rts(filename_prefix =file_prefix, savepath = f"{path_output}/",)
					obssch.write_ACPscript_LSGT(filename_prefix =file_prefix, period_script= 3, savepath = f"{path_output}/",) # for RASA36, shutdown = False
					logtbl = obssch.write_txt(filename_prefix =file_prefix, scheduled_only = False, format_ = 'ascii.fixed_width', savepath = f"{path_output}/",)
					obssch.show(save= False, filename_prefix = file_prefix, savepath = f"{path_output}/",)
					#------------------------------------------------------------
					#	LOAO
					#------------------------------------------------------------
					obssch = ObsScheduler(
						target_db=simple_galcat[:NUMBER_GALAXY_LIMIT],
						name_telescope='LOAO',
						entire_night=True,
						duplicate_when_empty=False,
						date=date,
						autofocus_at_init=False
						)
					obssch.write_rts(filename_prefix =file_prefix, savepath = f"{path_output}/",)
					logtbl = obssch.write_txt(filename_prefix =file_prefix, scheduled_only = False, format_ = 'ascii.fixed_width', savepath = f"{path_output}/",)
					obssch.show(save= False, filename_prefix = file_prefix, savepath = f"{path_output}/",)
					#------------------------------------------------------------
					#	SAO
					#------------------------------------------------------------
					obssch = ObsScheduler(
						target_db=simple_galcat[:NUMBER_GALAXY_LIMIT],
						name_telescope='SAO',
						entire_night=True,
						duplicate_when_empty=False,
						date=date,
						autofocus_at_init=False
						)
					obssch.write_rts(filename_prefix =file_prefix, savepath = f"{path_output}/",)
					logtbl = obssch.write_txt(filename_prefix =file_prefix, scheduled_only = False, format_ = 'ascii.fixed_width', savepath = f"{path_output}/",)
					obssch.show(save= False, filename_prefix = file_prefix, savepath = f"{path_output}/",)

				# %% [markdown]
				# # 4. sky grid matching

				# obs = 'KCT'
				# obs = 'KMTNet'
				# obs = 'CBNUO'
				# obs = 'RASA36'
				# for obs in ['KCT', 'CBNUO', 'RASA36']:
				for obs in ['KCT', 'CBNUO', 'RASA36', 'KMTNet',]:
					print(f"Generating Tiling Pattern for {obs}")

					pointing_cat = Table.read(f'../data/skygrid/{obs}/displaycenter.txt', format='csv')
					pointing_polygon_cat = Table.read(f'../data/skygrid/{obs}/displayfootprint.txt', format='csv')
					# coordinates = SkyCoord(pointing_cat['ra'], pointing_cat['dec'], unit='deg')

					# - count the number of pixels within 90% confidence area
					# 	- Full search 3.9s

					indx_select_skygrid = count_skymap_within_fov(pointing_polygon_cat, skymap, confidence_limit)
					select_pointing_cat = pointing_cat[indx_select_skygrid]
					select_pointing_polygon_cat = pointing_polygon_cat[indx_select_skygrid]

					from astropy.table import hstack
					select_skygrid_cat = hstack([select_pointing_cat, select_pointing_polygon_cat])
					if len(select_skygrid_cat) > 0:
						select_skygrid_cat['n_hostgalaxy']: int = 0
						select_skygrid_cat[probkey]: float = 0.0

						print(f"Number of points in selected pointing polygon: {len(select_pointing_polygon_cat)}")

						hostgalaxy90_point = PixCoord(simple_galcat['ra'], simple_galcat['dec'])

						for nn, (ra1, ra2, ra3, ra4, dec1, dec2, dec3, dec4) in enumerate(zip(select_skygrid_cat['ra1'], select_skygrid_cat['ra2'], select_skygrid_cat['ra3'], select_skygrid_cat['ra4'], select_skygrid_cat['dec1'], select_skygrid_cat['dec2'], select_skygrid_cat['dec3'], select_skygrid_cat['dec4'])):
							#	FoV Polygon
							vertices = PixCoord([ra1, ra2, ra3, ra4], [dec1, dec2, dec3, dec4],)
							region_pix = PolygonPixelRegion(vertices=vertices)

							n_hostgalaxy_points = len(simple_galcat[region_pix.contains(hostgalaxy90_point)])
							total_prob_polygon = np.sum(simple_galcat[probkey][region_pix.contains(hostgalaxy90_point)])
							# print(nn, n_hostgalaxy_points, total_prob_polygon)
							
							select_skygrid_cat['n_hostgalaxy'][nn] = n_hostgalaxy_points
							select_skygrid_cat[probkey][nn] = total_prob_polygon


						select_skygrid_cat.sort(probkey, reverse=True)
						for key in select_skygrid_cat.keys():
							if key in ['ra', 'dec', 'ra1', 'ra2', 'ra3', 'ra4', 'dec1', 'dec2', 'dec3', 'dec4', probkey]:
								select_skygrid_cat[key].format = '.3f'
						select_skygrid_cat['rank'] = np.arange(len(select_skygrid_cat))
						# select_skygrid_cat

						#   Cumulative Probability
						cumsum_prob_skygrid = np.cumsum(select_skygrid_cat[probkey])
						sum_prob_skygrid = np.sum(select_skygrid_cat[probkey])

						select_skygrid_cat['confidence'] = 0.0
						select_skygrid_cat['confidence'][np.max(cumsum_prob_skygrid)*1.0>=cumsum_prob_skygrid] = 1.0
						select_skygrid_cat['confidence'][np.max(cumsum_prob_skygrid)*0.95>=cumsum_prob_skygrid] = 0.95
						select_skygrid_cat['confidence'][np.max(cumsum_prob_skygrid)*0.9>cumsum_prob_skygrid] = 0.9
						select_skygrid_cat['confidence'][np.max(cumsum_prob_skygrid)*0.5>cumsum_prob_skygrid] = 0.5

						fig = plt.figure(figsize=(20, 8))
						fig = plt.figure(figsize=(40, 16))

						title = f"{record['alert_type']}: {area_90.to_value(u.deg**2):.1f} "+r"$\rm deg^2$ "+f"for {obs}"
						plt.subplot(121)
						if area_90.to(u.deg**2).value>50:
							only_center=True
						else:
							only_center=False
						plot_tiling_inorder(select_skygrid_cat, simple_galcat, skymap, title=title, only_center=only_center, probkey=probkey)
						xl, xr = plt.xlim()
						plt.xlim([xr, xl])
						plt.xticks(fontsize=20)
						plt.yticks(fontsize=20)
						plt.grid('both', color='silver', ls='--', lw=1, alpha=0.5)


						plt.subplot(122)
						# fig = plt.figure(figsize=(10, 8))
						plt.plot(cumsum_prob_skygrid, 'o-', mfc='w', mew=3, ms=10, lw=3, c='g')
						# include 95%, 99% 
						plt.axhline(y=0.95*sum_prob_skygrid, ls='-', c='tomato', lw=3, alpha=0.75, label=f"95% ({len(select_skygrid_cat[select_skygrid_cat['confidence']<=0.95])})")
						plt.axhline(y=0.9*sum_prob_skygrid, ls=':', c='tomato', lw=3, alpha=0.75, label=f"90% ({len(select_skygrid_cat[select_skygrid_cat['confidence']<=0.9])})")
						plt.axhline(y=0.5*sum_prob_skygrid, ls='--', c='tomato', lw=3, alpha=0.5, label=f"50% ({len(select_skygrid_cat[select_skygrid_cat['confidence']<=0.5])})")

						if probkey == 'prob_vol':
							ylabel = r'Cumulative $\rm P_{3D}$'
						elif probkey == 'prob_vol_x_stmass':
							ylabel = r'Cumulative $\rm P_{3D}xP_{M_{*}}$'
						plt.ylabel(ylabel)
						plt.xlabel("Number of Tiles")

						yl, yu = plt.ylim()
						plt.ylim([0, yu])

						plt.legend(loc='lower right', fontsize=20)
						# _ = plt.yticks(np.arange(0, 1.1, 0.1))
						plt.grid('both', color='silver', ls='--', lw=1, alpha=0.5)
						plt.title(f"{obs}")

						plt.grid('both', color='silver', ls='--', lw=1, alpha=0.5)
						plt.xticks(fontsize=20)
						plt.yticks(fontsize=20)


						plt.tight_layout()
						plt.savefig(f"{path_output}/tiling_{obs}.png", dpi=100)
						plt.close()

						cat_meta_dict = {
							"obs": obs,
							"alert_type": record['alert_type'],
							"superevent_id": record['superevent_id'],
						}
						select_skygrid_cat.meta = cat_meta_dict
						select_skygrid_cat_name = f"{path_output}/SkyGridCatalog_{obs}_{confidence_limit*1e2:g}.csv"
						select_skygrid_cat['obj'] = record['superevent_id']
						select_skygrid_cat['note'] = select_skygrid_cat['#id']
						select_skygrid_cat['weight'] = select_skygrid_cat[probkey]

						skygrid_meta_dict = {
							'superevent_id': record['superevent_id'],
							'alert_type': record['alert_type'],
							'most_probable_event': most_probable_event,
							'most_probable_event_prob': most_probable_event_prob,
							'confidence': confidence_limit,
							'ordering': probkey,
							"obs": obs,
						}
						select_skygrid_cat.write(f"{select_skygrid_cat_name}", format='csv', overwrite=True)

						date = Time.now()
						#	Observation Request
						if obs in ['RASA36', 'KCT', 'CBNUO']:
							data_original = select_skygrid_cat
							# data_original.remove_rows([7, 4])
							subsequent_days = obs_request_day
							file_prefix ='GECKO_'
							for i in range(subsequent_days):
								obssch = ObsScheduler(target_db = data_original, name_telescope = obs, entire_night = True, duplicate_when_empty= False, date = date, autofocus_at_init= False)
								if obs == 'RASA36':
									obssch.write_ACPscript_RASA36(filename_prefix =file_prefix, savepath = f"{path_output}/",) # for RASA36, shutdown = False   
								elif obs == 'KCT':
									obssch.write_ACPscript_KCT(filename_prefix =file_prefix, shutdown = True, savepath = f"{path_output}/",) # for RASA36, shutdown = False
								else:
									obssch.write_rts(filename_prefix =file_prefix, savepath = f"{path_output}/")
								# elif obs == 'LSGT':
								# 	obssch.write_ACPscript_LSGT(filename_prefix =file_prefix, period_script= 3) # for RASA36, shutdown = False
								logtbl = obssch.write_txt(filename_prefix =file_prefix, scheduled_only = False, format_ = 'ascii.fixed_width', savepath = f"{path_output}/",)
								obssch.show(save= False, filename_prefix = file_prefix, savepath = f"{path_output}/",)
								plt.close()
								date += 1*u.day
						elif obs == 'KMTNet':
							data_original = select_skygrid_cat
							for i in range(subsequent_days):
								cbands = ['B', 'R']
								date_str = f"{date.datetime.year}{date.datetime.month:0>2}{date.datetime.day:0>2}"

								for site in ['SSO', 'SAAO', 'CTIO']:
									obsname = f"KMTNet_{site}"
									file_prefix ='GECKO_'
									obssch = ObsScheduler(target_db = data_original, name_telescope = obsname, entire_night = True, duplicate_when_empty= False, date = date, autofocus_at_init= False)
									logtbl = obssch.write_txt(filename_prefix =file_prefix, scheduled_only = False, format_ = 'ascii.fixed_width', savepath = f"{path_output}/",)
									obssch.show(save= False, filename_prefix = file_prefix, savepath = f"{path_output}/",)
									plt.close()

									#	KMTNet Script Inputs
									cfield = [record['superevent_id']]*len(logtbl)
									crarr = logtbl['ra']
									cdecarr = logtbl['dec']

									generate_KMTNet_obs_script(
										cfield,
										crarr,
										cdecarr,
										site,
										date_str,
										cbands,
										path_save=path_output)
								date += 1*u.day

				# %% [markdown]
				# # Fin. ALERT!!!

				# %%
				delt = time.time() - st
				OAuth_Token = SLACK_API_CONFIG['OAuth_Token']

				channel = '#gecko--alert'
				text = f"[`GeckoDigestor`] Observation is ready for {record['superevent_id']}-{record['alert_type']} ({delt:1.1f}s)!\nqso:{path_output}"
				#---------------------------------
				#       Slack message
				#---------------------------------
				param_slack = dict(
						token = OAuth_Token,
						channel = channel,
						text = text,
				)
				if slack:
					slack_bot(**param_slack)
			else:

				OAuth_Token = SLACK_API_CONFIG['OAuth_Token']

				channel = '#gecko--alert'
				text = f"[`GeckoDigestor`] Catch the alert ({record['superevent_id']}-{record['alert_type']}). Stay tuned!"
				#---------------------------------
				#       Slack message
				#---------------------------------
				param_slack = dict(
						token = OAuth_Token,
						channel = channel,
						text = text,
				)
				if slack:
					slack_bot(**param_slack)
		else:
			OAuth_Token = SLACK_API_CONFIG['OAuth_Token']

			channel = '#gecko--alert'
			text = f"[`GeckoDigestor`] {record['superevent_id']}-{record['alert_type']}"
			#---------------------------------
			#       Slack message
			#---------------------------------
			param_slack = dict(
					token = OAuth_Token,
					channel = channel,
					text = text,
			)
			if slack:
				slack_bot(**param_slack)
