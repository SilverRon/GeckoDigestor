from astropy.coordinates import SkyCoord
from regions import PixCoord, PolygonSkyRegion, PolygonPixelRegion
from regions import RegularPolygonPixelRegion
import numpy as np

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