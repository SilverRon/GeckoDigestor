from astropy.coordinates import SkyCoord
from regions import PixCoord, PolygonSkyRegion, PolygonPixelRegion
from regions import RegularPolygonPixelRegion
import numpy as np
import healpy as hp
from astropy_healpix import uniq_to_level_ipix

def count_skymap_within_fov(skygrid_cat, skymap, confidence_limit):
	# count_skymap_within_fov(pointing_polygon_cat, skymap, confidence_limit)
	#	90% --> 20.2s
	skymap90 = skymap[skymap['CUMPROBDENSITY']<confidence_limit]
	#	100% --> 40.7s
	# skymap90_points = PixCoord(skymap90['RA'], skymap90['DEC'])

	indx_select_skygrid = []
	for nn, (ra, ra1, ra2, ra3, ra4, dec, dec1, dec2, dec3, dec4) in enumerate(zip(skygrid_cat['ra'], skygrid_cat['ra1'], skygrid_cat['ra2'], skygrid_cat['ra3'], skygrid_cat['ra4'], skygrid_cat['dec'], skygrid_cat['dec1'], skygrid_cat['dec2'], skygrid_cat['dec3'], skygrid_cat['dec4'])):
		if ra == 0.0:
			ra3 -= 360
			ra4 -= 360
		else:
			pass
		radistance = max([ra1, ra2, ra3, ra4,])-min([ra1, ra2, ra3, ra4,])
		ramin = min([ra1, ra2, ra3, ra4])
		ramax = max([ra1, ra2, ra3, ra4])
		decmin = min([dec1, dec2, dec3, dec4])
		decmax = max([dec1, dec2, dec3, dec4])
		_skymap90 = skymap90[
			(skymap90['RA']>=ramin-5) &
			(skymap90['RA']<=ramax+5) &
			(skymap90['DEC']>=decmin) &
			(skymap90['DEC']<=decmax)
			].copy()
		# if (radistance < 100):
		# if (len(_skymap90) > 0) & (radistance < 150):
		if (len(_skymap90) > 0):
			skymap90_points = PixCoord(_skymap90['RA'], _skymap90['DEC'])

			vertices = PixCoord([ra1, ra2, ra3, ra4], [dec1, dec2, dec3, dec4],)
			region_pix = PolygonPixelRegion(vertices=vertices)
			number_of_skymap_points = len(_skymap90[region_pix.contains(skymap90_points)])

			if (number_of_skymap_points > 0):
				indx_select_skygrid.append(nn)
	#	List --> numpy array
	indx_select_skygrid = np.array(indx_select_skygrid)
	return indx_select_skygrid

# def count_skymap_within_fov_fast(skygrid_cat, skymap, confidence_limit, nside):
#     # 1) 90% 영역 픽셀 UNIQ 세트
#     keep = skymap['CUMPROBDENSITY'] < confidence_limit
#     uniq90 = set(skymap['UNIQ'][keep])

#     sel_idx = []
#     for i, row in enumerate(skygrid_cat):
#         # 2) 타일 사각형의 (lon,lat) → 벡터
#         lons = np.radians([row['ra1'],row['ra2'],row['ra3'],row['ra4']])
#         lats = np.radians([row['dec1'],row['dec2'],row['dec3'],row['dec4']])
#         vec  = hp.ang2vec(0.5*np.pi - lats, lons)

#         # 3) 폴리곤 내 모든 픽셀 UNIQ 리스트
#         pixs = hp.query_polygon(nside, vec, inclusive=True)

#         # 4) 교집합 검사
#         if uniq90.intersection(pixs):
#             sel_idx.append(i)

#     return np.array(sel_idx)

def count_skymap_within_fov_fast(skygrid_cat, skymap, confidence_limit):
    """
    빠른 버전: astropy_healpix(ah)와 healpy(hp)를 조합해서
    90% 영역을 포함하는 타일 인덱스만 골라냅니다.
    """
    # 1) skymap['UNIQ'] → (level, ipix_map) & nside_map
    level, ipix_map = ah.uniq_to_level_ipix(skymap['UNIQ'])
    nside_map        = ah.level_to_nside(level)   # 보통 전체가 같은 레벨이라 하나만 취해도 무방

    # 2) 90% 영역에 속한 ipix 집합
    mask90 = skymap['CUMPROBDENSITY'] < confidence_limit
    ipix90 = set(ipix_map[mask90])

    sel_idx = []
    for i, row in enumerate(skygrid_cat):
        # 3) 타일 꼭짓점 → 천구 벡터
        lons = np.radians([row['ra1'], row['ra2'], row['ra3'], row['ra4']])
        lats = np.radians([row['dec1'], row['dec2'], row['dec3'], row['dec4']])
        vecs = hp.ang2vec(0.5*np.pi - lats, lons)

        # 4) 같은 nside_map 해상도로 타일 영역 픽셀 ipix 뽑기
        tile_ipix = hp.query_polygon(nside_map, vecs, inclusive=True)

        # 5) 90% 영역 ipix와 교집합이 있으면 ‘관통’
        if ipix90.intersection(tile_ipix):
            sel_idx.append(i)

    return np.array(sel_idx, dtype=int)


def is_condition_satisfied(record):
    is_classification_NS = 'NS' in record['event']['classification']
    is_burst = 'Burst' in record['event']['group']
    is_distmean_less = record['distmean'] <= 1_000 # [Mpc]
    is_area_NS_less = record['area_90'] <= 15_000 # [deg2]
    is_area_BH_less = record['area_90'] <= 5_000 # [deg2]
    is_area_BURST_less = record['area_90'] <= 100 # [deg2]

    return (is_classification_NS and is_distmean_less and is_area_NS_less) or (
        not is_classification_NS and is_distmean_less and is_area_BH_less
    ) or (
	    is_burst and is_distmean_less and is_area_BURST_less
	)

def DrawTiles(intbl):
	allralist = []
	alldeclist = []
	for jj in np.arange(len(intbl)):
		ralist = []
		declist = []
		for ii in [1, 2, 3, 4, 1]:
			ra = intbl[f'ra{ii}'][jj]
			dec = intbl[f'dec{ii}'][jj]
			ralist.append(ra)
			declist.append(dec)
		allralist.append(ralist)
		alldeclist.append(declist)
	return allralist, alldeclist


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

from astropy.coordinates import SkyCoord
import astropy.units as u

def makeSpecColors(n, palette='Spectral'):
	#	Color palette
	palette = sns.color_palette(palette, as_cmap=True,)
	palette.reversed

	clist_ = [palette(i) for i in range(palette.N)]
	cstep = int(len(clist_)/n)
	clist = [clist_[i*cstep] for i in range(n)]
	return clist


def plot_tiling_inorder(select_skygrid_cat, simple_galcat, skymap, title="", only_center=False, probkey='prob_vol'):
	#	Beautiful colors
	colors = makeSpecColors(len(select_skygrid_cat), palette="cividis")

	#
	cumprob = skymap['CUMPROBDENSITY']

	# fig = plt.figure(figsize=(10, 8))
	#	Galaxy
	if only_center:
		markersize = 6
	else:
		markersize = 36
	if not only_center:
		indx_sort = np.argsort(simple_galcat[probkey])
		plt.scatter(simple_galcat['ra'][indx_sort], simple_galcat['dec'][indx_sort], c=simple_galcat[probkey][indx_sort], cmap='hot', edgecolors='k', s=markersize, alpha=0.75, zorder=10)
		cbar_gal = plt.colorbar()

	plt.scatter(skymap['RA'][cumprob<0.9], skymap['DEC'][cumprob<0.9], c=skymap['CUMPROBDENSITY'][cumprob<0.9], alpha=0.5)
	# cbar = plt.colorbar()

	# plt.xlim([xr, xl])
	# plt.ylim([yl, yu])
	if only_center:
		plt.plot(select_skygrid_cat['ra'], select_skygrid_cat['dec'], '+', c='red')

	for nn in np.arange(len(select_skygrid_cat)):
		ralist = [0.0]*5
		declist = [0.0]*5
		if only_center == False:
			plt.text(select_skygrid_cat['ra'][nn], select_skygrid_cat['dec'][nn], nn, ha='center', va='center', fontsize=14, zorder=11)
		for ii, jj in enumerate([1, 2, 3, 4, 1]):
			ra, dec = select_skygrid_cat[f'ra{jj}'][nn], select_skygrid_cat[f'dec{jj}'][nn]
			ralist[ii] = ra
			declist[ii] = dec
		color = colors[nn]
		if nn != 0:
			plt.plot(ralist, declist, '-', c=color)
		else:
			plt.plot(ralist, declist, '-', c=color, label=f'Tiles: {len(select_skygrid_cat)}')

	xl, xr = plt.xlim()
	yl, yu = plt.ylim()

	plt.xlim([xl, xr])
	plt.ylim([yl, yu])

	# plt.title(f"{record['alert_type']}: {area_90.to_value(u.deg**2):.1f} "+r"$\rm deg^2$")
	plt.title(title)
	plt.xlabel('RA [deg]')
	plt.ylabel('Dec [deg]')
	plt.legend(fontsize=14)


def min_max_normalize(data):
    """
    Normalize a given dataset using min-max normalization
    """
    min_val = min(data)
    max_val = max(data)
    normalized_data = [(x - min_val) / (max_val - min_val) for x in data]
    return np.array(normalized_data)


import numpy as np
from astropy import units as u
import astropy.coordinates as coord

def generate_KMTNet_obs_script(cfield: list, craarr: np.ndarray, cdecarr: np.ndarray, site: str, date: str, cbands: list, cid: list, path_save: str, cdith=[0, 3]):
	'''
	This function generates a KMTNet observation script given inputs for object name, RA/Dec coordinates, observation site, date, and filter bands.
	The function takes the following arguments:

	Arguments:
	- `cfield`: An array of object names.
	- `craarr`: An array of right ascension (RA) values in degrees or in string format (e.g., '19:08:00.000') corresponding to each object in `cfield`.
	- `cdecarr`: An array of declination (Dec) values in degrees or in string format (e.g., '-64:30:00.00') corresponding to each object in `cfield`.
	- `site`: The name of the observation site. It can be 'SSO', 'SAAO', or 'CTIO'.
	- `date`: The date of the observation in the format 'YYYYMMDD'.
	- `cbands`: An array of filter bands to be used for the observation.
	
	The function returns nothing, but creates a file named 'scTOO_YYYYMMDD_site.cat', which contains the observation script in a specific format.
	The observation script includes the project ID, label, RA and Dec coordinates, imaging type, object name, filter, exposure time, UTC start time, and UTC tolerance for each object and filter combination.
	The script is created by iterating through each object in `cfield`, each dithering option (0=no dithering, 1=RA direction, 2=Dec direction, 3=RA&Dec direction), and each filter in `cbands`.
	The function includes inline comments to explain each step of the observation script generation process.
	'''
	#------------------------------------------------------------
	#	KMTNet Observation Script
	#------------------------------------------------------------
	#	2023.04.21 Modulizied by Gregory S.H. Paek
	#	202?.??.?? Generated by Seo-won Chang
	#------------------------------------------------------------
	#cfield = object name array
	#craarr = ra[deg] array
	#cdecarr = dec[deg] array
	#cband = filter array
	#cdith = dithering array 0,1,2,3 (0=no dithering, 1=ra direction, 2=dec direction, 3=ra&dec direction)
	#------------------------------------------------------------#------------------------------------------------------------
	# 2022B DWF
	# NGC6744         19:08:00.000  -64:30:00.00          04:45 - 06:30 # For NGC 6744, the field is offset from the center of the galaxy
	# FRB190711         21:57:40.680  -80:21:28.80          06:30 - 09:00 # neeed to use off-set coordinates from pts either as before
	# HDF-S                   22:33:26.000  -60:38:09.00          09:00 - 09:30
	# CDF-S                    03:30:00.000  -28:06:00.00         09:30 - 10:15

	# output example
	# KMTNet Observation script
	#  ProjectID     LABEL         RA          DEC     COPT  IMGTYP   OBJECT_NAME FILTER EXPTIME      UTOBS          UTTOL #Comments
	#------------ ----------- ----------- ---- -------- ------------ ------ ------- ------------------- ----- ---------
	# TOO         FRB171216-0(V) 03:28:00.0   -57:04:00    0   OBJECT    FRB171216-0    V     120                   -     0 #1
	#------------------------------------------------------------#------------------------------------------------------------
	# site = 'SSO'
	# site = 'SAAO'
	# site = 'CTIO'

	# site = 'CTIO'
	# date = '20230313' # 16
	# cdith = [0,1,2,3, 4,5,6,7, 8,9,10,11, 12,13,14,15, 16,17,18,19]
	# cdith = [0, 3]
	# cbands = ['R',]
	#craarr = [287.000000, 329.419500, 338.358333, 52.500000]
	#cdecarr = [-64.500000, -80.358000, -60.635833, -28.100000]
	#craarr = ['19:08:00.000', '21:57:40.680', '22:33:26.000', '03:30:00.000']
	#cdecarr = ['-64:30:00.00', '-80:21:28.80', '-60:38:09.00', '-28:06:00.00']

	#------------------------------------------------------------
	#	Input Data
	# craarr = [61.277]
	# cdecarr = [-76.000]
	# cfield = ['TOO_0503',]
	#------------------------------------------------------------
	f = open(f'{path_save}/scTOO_'+date+'_'+site+'.cat',"w")
	f.write('# KMTNet Observation script\n')
	f.write('#  ProjectID     LABEL         RA          DEC     COPT  IMGTYP   OBJECT_NAME FILTER EXPTIME      UTOBS          UTTOL #Comments\n')
	f.write('#------------ ----------- ----------- ---- -------- ------------ ------ ------- ------------------- ----- ---------\n')
	f.write('\n')

	cnt = 0
	for i in range(len(craarr)):
		# iteration by field name
		tname = cfield[i]
		cidname = cid[i]
		
		# iteration by dithing option
		# (0=no dithering, 1=ra direction, 2=dec direction, 3=ra&dec direction)
		for dith in cdith:
			tdra = craarr[i]
			tddec = cdecarr[i]
			if dith == 1: 
				tdra = craarr[i] + ((4/60.) / np.cos(cdecarr[i] * np.pi / 180.))
			if dith == 2:
				tddec = cdecarr[i] + 7/60.
			if dith == 3: 
				tdra = craarr[i] + ((4/60.) / np.cos(cdecarr[i] * np.pi / 180.))
				tddec = cdecarr[i] + 7/60.

			tra_ = coord.Angle(tdra, unit=u.deg)
			tra_sep = tra_.to_string(unit=u.hour, sep=':',precision=1, pad=True)
			tdec_ = coord.Angle(tddec, unit=u.deg)
			tdec_sep = tdec_.to_string(unit=u.degree, sep=':',precision=0, pad=True)

			# iteration by filters
			for cband in cbands:
				cnt = cnt + 1
				f.write('{0:<12}{1:<15}{2:<13}{3:<13}{4:<14}{5:<15}{6:<6}{7:<29}{8:}\n'.format('TOO',tname+'-'+str(dith), tra_sep,tdec_sep,'0   OBJECT',tname,cband,'120                   -     0',' #'+str(cidname)))
	f.close()
	return 


import numpy as np

def calc_apparent_mag(d, dprime, m):
    """
    Calculate the apparent magnitude at a distance dprime, given the apparent magnitude m at distance d.

    Parameters:
    -----------
    d : float
        Distance to the object in parsecs
    dprime : float
        Distance to the object in parsecs for which the magnitude is to be calculated
    m : float
        Apparent magnitude of the object at distance d

    Returns:
    --------
    mprime : float
        Apparent magnitude of the object at distance dprime
    """
    mprime = m + 5 * np.log10(dprime / d)
    return mprime


def expect_AT2017gfo(dprime, phase, path_template_lc, filterlist=['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K',], plot=True):
	# filterlist = ['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K',]
	# for filte in ['g', 'r', 'i', 'z', 'y',]:
	colors = makeSpecColors(len(filterlist), palette='Spectral')
	expected_magdict = {}
	d = 40 # [Mpc]
	#	AT2017gfo Table
	# intbl = ascii.read('../data/AT2017gfo_phot_modified.dat', header_start=0, data_start=1)
	intbl = ascii.read(path_template_lc, header_start=0, data_start=1)
	for ff, filte in enumerate(filterlist):
		# plt.plot(intbl['Phase'], intbl[filte], 'o-', mfc='w', label=filte)
		mprime = calc_apparent_mag(d, dprime, intbl[filte])
		mag_now = np.interp(phase, intbl['Phase'][~intbl[filte].mask], mprime[~intbl[filte].mask])
		# if mag_now > 10:
		# 	label = f"{filte}={mag_now:.1f}"
		# else:
		# 	label = filte
		label = f"{filte}={mag_now:.1f}"
		expected_magdict[filte] = mag_now
		if plot:
			# plt.plot(intbl['Phase'], mprime, 'o-', mfc='w', label=label)
			# plt.plot(intbl['Phase'][~mprime.mask], mprime[~mprime.mask], 'o-', mfc='w', label=label)
			plt.plot(intbl['Phase'][~mprime.mask], mprime[~mprime.mask], 'o-', c=colors[-ff-1], mfc='w', label=label)

	if plot:
		plt.axvspan(xmin=1e-3, xmax=phase, color='grey', alpha=0.25)
		plt.axvline(x=phase, ls='-', color='tomato', lw=3, alpha=1.0, zorder=1, label=f'Now {phase:.1f}d')
		plt.xlabel('Phase')
		plt.ylabel('Mag [AB]')
		plt.title(f"AT2017gfo-like KN (d={dprime:.1f} Mpc)")
		plt.xlim([0, phase+2])
		# yl, yu = plt.ylim()
		yl = np.min(list(expected_magdict.values()))-0.25
		yu = np.max(list(expected_magdict.values()))+3
		plt.ylim([yu, yl])
		plt.legend(framealpha=0.0, fontsize=10, loc='lower left', ncol=3)
		plt.grid('both', ls='--', c='silver', alpha=0.5)
		plt.tight_layout()
	return expected_magdict


def scale_depth(depth, t0, t1):
    tratio = t1/t0
    return depth+2.5*np.log10(np.sqrt(tratio))

def scale_exptime(m0, m1):
    delm = m1 - m0
    tratio = 10**(delm/(1.25))
    return tratio



def deg2hmsdms(ra_deg, dec_deg):
    """
    Convert RA, Dec in decimal degrees to sexagesimal HMS, DMS strings.
    
    Parameters
    ----------
    ra_deg : float or numpy.ndarray
        Right Ascension in decimal degrees
    dec_deg : float or numpy.ndarray
        Declination in decimal degrees
    
    Returns
    -------
    hms, dms : tuple of str or tuple of numpy.ndarray
        Sexagesimal HMS string and DMS string
    """
    # Create SkyCoord object with input RA, Dec in degrees
    coords = SkyCoord(ra_deg*u.deg, dec_deg*u.deg, frame='icrs')
    # Convert to sexagesimal format
    hms = coords.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
    dms = coords.dec.to_string(unit=u.degree, sep=':', precision=2, pad=True, alwayssign=True)
    
    return hms, dms


def exptime_for_mag(m, depth, exptime0):
    """
    Given the brightness of the target source m, the depth of the telescope, and the initial exposure time,
    returns the exposure time required to observe the target source.

    Parameters
    ----------
    m : float or numpy.ndarray
        The brightness of the target source in magnitudes.
    depth : float
        The depth of the telescope in magnitudes.
    exptime0 : float
        The initial exposure time in seconds.

    Returns
    -------
    float or numpy.ndarray
        The exposure time required to observe the target source in seconds.
    """
    exptime = exptime0 * 10**((m - depth) / 2.5)
    return exptime

def find_exposure_time(exp_min, obs_min, total_time, tolerence=10):
	closest_time = float('inf')
	closest_tuple = None

	for n in range(obs_min, int(total_time // exp_min) + 1):
		for t in range(exp_min, 301, 30):
			exp_time = n * t
			diff = abs(exp_time - total_time)
			if diff < closest_time:
				closest_time = diff
				closest_tuple = (t, n)
	if closest_tuple == None:
		closest_tuple = (exp_min, obs_min)
	return closest_tuple

def create_ds9_regions(ra, dec, names, filename='regions_with_names.reg'):
	"""
	Create a DS9 region file with circles of radius 5 arcseconds and include object names as labels.
	Correctly formats text labels to ensure compatibility with DS9.

	Parameters:
	- ra: array of Right Ascension in degrees
	- dec: array of Declination in degrees
	- names: array of object names
	- filename: name of the output region file
	"""
	with open(filename, 'w') as file:
		# Write the DS9 region file header
		file.write('# Region file format: DS9 version 4.1\n')
		file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" ')
		file.write('select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
		file.write('fk5\n')

		# Write each circle region with the object name
		for r, d, name in zip(ra, dec, names):
			file.write(f'circle({r},{d},5") # text={{{name}}}\n')

	print(f"DS9 region file '{filename}' created successfully.")