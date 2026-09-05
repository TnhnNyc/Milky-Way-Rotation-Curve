from __future__ import annotations
import numpy
from astropy.table import join

# Adjust these if your downloaded APOGEE data set uses different column names
COLUMN_MAP = {
    "gaia_source_id": "source_id",
    "distance_pc": "weighted_dist",
    "distance_pc_err": "weighted_dist_error",
    "logg": "LOGG",
    "teff": "TEFF",
    "feh": "FE_H",
    "ra": "RA",
    "dec": "DEC"
    }

#===============================================
# Filtering Dataset according to the Ruwe Value
#===============================================
##(See the Milky_Way_Mass_RotationCurve.ipynb notebook, Section 1.2)

def quality_cuts(dataset, ruwe_max, parallax_snr=None):
    """
    Basic Gaia astrometric quality cut. RUWE < 1.4 is the standard rule
    of thumb for a well-behaved single-star astrometric solution (Lindegren
    et al. 2021). Since we're using astroNN spectrophotometric distances
    rather than Gaia parallax for distance, a parallax S/N cut is optional
    here -- but keeping RUWE clean still matters for proper motions.
    """
    mask = (
        (dataset['ruwe'] < ruwe_max)
        & (dataset[COLUMN_MAP["distance_pc"]] > 0)
    )

    if parallax_snr is not None:
        snr = dataset['nn_parallax'] / dataset['nn_parallax_error']
        mask &= snr > parallax_snr
        
    return dataset[mask]

#===============================================
# Selecting the Red Giant Stars
#===============================================
##(See the Milky_Way_Mass_RotationCurve.ipynb notebook, Section 1.2)

def select_red_giants(dataset, logg_max, teff_min, teff_max):
    mask = (
        (dataset['LOGG'] < logg_max)
        & (dataset['TEFF'] < teff_max)
        & (dataset['TEFF'] > teff_min)
    )
    return dataset[mask]

#===============================================
# Selecting the Outer Disk Red Giant Stars
#===============================================
##(See the Milky_Way_Mass_RotationCurve.ipynb notebook, Section 2.2)

def select_outer_disk(dataset, R_kpc, z_kpc, R_min, R_max, z_max):
    """
    Apply the R > 15 kpc outer-disk cut (and an optional |z| cut to
    reduce halo contamination -- outer-disk *disk* stars should stay
    reasonably close to the plane even though the disk flares at large R).

    R_kpc, z_kpc are the arrays returned by
    kinematics.galactocentric_cylindrical_velocities -- compute those
    first, then call this to filter both the table and the kinematics.
    """
    mask = (R_kpc >= R_min) & (R_kpc <= R_max) & (numpy.abs(z_kpc) <= z_max)
    return mask
