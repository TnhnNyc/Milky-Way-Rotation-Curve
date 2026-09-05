from __future__ import annotations
import numpy
import astropy.units as units
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, Galactocentric



# ======================================================================
# Calculation of the Galactocentric 3D Velocities of the Red Giant Stars 
# ======================================================================

# 1a. Creating a galactocentric frama: 
def get_galactocentric_frame(galcen_distance=None, z_sun=None, galcen_v_sun=None):
    """
    Build a Galactocentric frame. With no arguments this uses astropy's
    current IAU-recommended defaults. Pass explicit values if you want
    to match a specific paper's assumed solar position/velocity exactly
    (e.g. to reproduce Jiao et al. 2023's numbers).
    """
    kwargs = {}
    if galcen_distance is not None:
        kwargs["galactocentric_distance"] = galcen_distance
    if z_sun is not None:
        kwargs["z_sun"] = z_sun
    if galcen_v_sun is not None:
        kwargs["galactocentric_v_sun"] = galactocentric_v_sun
    frame = Galactocentric(**kwargs)
    return frame

# 1b. Calculating the 3D velocities: 
def observables_to_galactocentric(ra, dec, distance, pmra, pmdec, radial_velocity, frame=None):
    """
    Parameters
    ----------
    ra, dec          : Quantity [deg]
    distance         : Quantity [kpc]
    pmra, pmdec      : Quantity [mas/yr]
    radial_velocity  : Quantity [km/s]
    frame            : astropy Galactocentric frame (optional)
    """
    if frame is None:
        frame = get_galactocentric_frame()

    icrs = SkyCoord(
        ra=ra, dec=dec, distance=distance,
        pm_ra_cosdec=pmra, pm_dec=pmdec,
        radial_velocity=radial_velocity,
        frame='icrs'
    )
    return icrs.transform_to(frame)

# =====================================================================================
# Calculation the Velocities of the Red Giant Stars in Cylindrical Coordinate System
# =====================================================================================

def galactocentric_cylindrical_velocities(cartesian_velocities):
    """
    Convert Galactocentric Cartesian coordinates and velocities
    into cylindrical coordinates and velocity components.
    """

    # Cartesian positions and velocities
    x = cartesian_velocities.x.to_value(units.kpc)
    y = cartesian_velocities.y.to_value(units.kpc)
    z = cartesian_velocities.z.to_value(units.kpc)
    vx = cartesian_velocities.v_x.to_value(units.km / units.s)
    vy = cartesian_velocities.v_y.to_value(units.km / units.s)
    vz = cartesian_velocities.v_z.to_value(units.km / units.s)

    #Cylindrical position
    R = numpy.sqrt(x**2 + y**2)
    phi = numpy.arctan2(y,x)

    #Cylindrical velocity components
    v_R = (x * vx + y * vy) / R
    v_phi_math = (x * vy - y * vx) / R

    #Astropy Galactocentric convention
    v_phi = -v_phi_math

    return {
        "R_kpc":R,
        "phi_rad":phi,
        "z_kpc":z,
        "v_R_kms":v_R,
        "v_phi_kms":v_phi,
        "v_z_kms":vz
    }
    
