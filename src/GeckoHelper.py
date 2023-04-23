from astropy.coordinates import SkyCoord
from regions import PixCoord, PolygonSkyRegion, PolygonPixelRegion
from regions import RegularPolygonPixelRegion
import numpy as np

def count_skymap_within_fov(pointing_polygon_cat, skymap, confidence_limit):
	# count_skymap_within_fov(pointing_polygon_cat, skymap, confidence_limit)
	#	90% --> 20.2s
	skymap90 = skymap[skymap['CUMPROBDENSITY']<confidence_limit]
	#	100% --> 40.7s
	skymap90_points = PixCoord(skymap90['RA'], skymap90['DEC'])

	indx_select_skygrid = []
	for nn, (ra1, ra2, ra3, ra4, dec1, dec2, dec3, dec4) in enumerate(zip(pointing_polygon_cat['ra1'], pointing_polygon_cat['ra2'], pointing_polygon_cat['ra3'], pointing_polygon_cat['ra4'], pointing_polygon_cat['dec1'], pointing_polygon_cat['dec2'], pointing_polygon_cat['dec3'], pointing_polygon_cat['dec4'])):
		radistance = max([ra1, ra2, ra3, ra4,])-min([ra1, ra2, ra3, ra4,])

		if (radistance < 100):
			vertices = PixCoord([ra1, ra2, ra3, ra4], [dec1, dec2, dec3, dec4],)
			region_pix = PolygonPixelRegion(vertices=vertices)
			number_of_skymap_points = len(skymap90[region_pix.contains(skymap90_points)])

			if (number_of_skymap_points > 0):
				indx_select_skygrid.append(nn)
	#	List --> numpy array
	indx_select_skygrid = np.array(indx_select_skygrid)
	return indx_select_skygrid