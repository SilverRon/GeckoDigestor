#%%
#   Library
import astropy.coordinates
import astropy.time
import astropy.units as u
import healpy as hp
#%%
#   Function
def prob_observable(m, header):
    """
    Determine the integrated probability contained in a gravitational-wave
    sky map that is observable from a particular ground-based site at a
    particular time.

    Bonus: make a plot of probability versus UTC time!
    """

    # Determine resolution of sky map
    npix = len(m)
    nside = hp.npix2nside(npix)

    # Get time now
    time = astropy.time.Time.now()
    # Or at the time of the gravitational-wave event...
    # time = astropy.time.Time(header['MJD-OBS'], format='mjd')
    # Or at a particular time...
    # time = astropy.time.Time('2015-03-01 13:55:27')

    # Geodetic coordinates of observatory (example here: Mount Wilson)
    observatory = astropy.coordinates.EarthLocation(
        lat=34.2247*u.deg, lon=-118.0572*u.deg, height=1742*u.m)

    # Alt/az reference frame at observatory, now
    frame = astropy.coordinates.AltAz(obstime=time, location=observatory)

    # Look up (celestial) spherical polar coordinates of HEALPix grid.
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    # Convert to RA, Dec.
    radecs = astropy.coordinates.SkyCoord(
        ra=phi*u.rad, dec=(0.5*np.pi - theta)*u.rad)

    # Transform grid to alt/az coordinates at observatory, now
    altaz = radecs.transform_to(frame)

    # Where is the sun, now?
    sun_altaz = astropy.coordinates.get_sun(time).transform_to(altaz)

    # How likely is it that the (true, unknown) location of the source
    # is within the area that is visible, now? Demand that sun is at
    # least 18 degrees below the horizon and that the airmass
    # (secant of zenith angle approximation) is at most 2.5.
    prob = m[(sun_altaz.alt <= -18*u.deg) & (altaz.secz <= 2.5)].sum()

    # Done!
    return prob
