from __future__ import annotations
import numpy
from astropy.table import join, Column
import astropy.units as units

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
# Removing Duplicate Star Entries
#===============================================
##(See the Milky_Way_Mass_RotationCurve.ipynb notebook, Section 1.1)

def deduplicate_stars(dataset, id_column="APOGEE_ID", priority_column="TEFF_ERR"):
    """
    A star can be appeared more than one under the same APOGEE_ID but different
    LOCATION_ID. It causes to increase star count and over-weight of multiply-
    observed stars. This function keeps the "best" row per unique star.

    The best row is identified by the lowest TEFF_ERR value. Invalid priority
    values (NaN or <= 0) are set to +inf to exclude corrupt data.
    """
    ids = numpy.asarray(dataset[id_column])
    priority = numpy.asarray(dataset[priority_column], dtype=float)
    priority = numpy.where(numpy.isnan(priority) | (priority <= 0), numpy.inf, priority)

    # Sort all rows by priority (best/lowest first); for each unique ID,
    # its first appearance in that ordering is automatically its best duplicate.
    order = numpy.argsort(priority, kind="stable")
    _, first_occurrence = numpy.unique(ids[order], return_index=True)
    keep_idx = numpy.sort(order[first_occurrence])

    n_before = len(dataset)
    result = dataset[keep_idx]
    n_after = len(result)
    print(
        f"Deduplication by {id_column}: {n_before} rows -> {n_after} unique "
        f"stars ({n_before - n_after} duplicate entries removed)"
    )
    return result

#===============================================
# Filtering Dataset according to the Ruwe Value
#===============================================
##(See the Milky_Way_Mass_RotationCurve.ipynb notebook, Section 1.2)

def quality_cuts(dataset, ruwe_max, parallax_snr=None, max_dist_err=None):
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

    if max_dist_err is not None:
        distance = dataset[COLUMN_MAP["distance_pc"]]
        distance_error = dataset[COLUMN_MAP["distance_pc_err"]]
        rel_error = distance_error / distance
        mask &= rel_error < max_dist_err 

    n_before = len(dataset)
    result = dataset[mask]
    print(
        f"quality_cuts: {n_before} -> {len(result)} stars "
        f"(RUWE<{ruwe_max}"
        + (f", dist_rel_err<{max_dist_err }" if max_dist_err  is not None else "")
        + (f", parallax_snr>{parallax_snr}" if parallax_snr is not None else "")
        + ")"
    )
    
    return result

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
# Filling Missing Gaia Radial Velocities
#===============================================
##(See the Milky_Way_Mass_RotationCurve.ipynb notebook, Section 1.2)

def fill_missing_radial_velocity(dataset,
                                 gaia_column,
                                 gaia_error_column,
                                 apogee_column,
                                 apogee_error_column=None,
                                 apogee_manual_error=0.5,
                                 sentinel_abs_limit=1000.0):
    """
    Gaia DR3 can fall behind for faint/distant giants, and it leads to
    missing data (Nan) in a data table. This causes some issues:
    the rows including Nan are not put into calculation of velocity
    moments, but they are treated as one star in the star count.

    To resolve this, APOGEE's spectroscopic RV (shown as`VHELIO_AVG`
    in the data), is used as a fallback, since it covers radial
    velocity data  for almost all giant.

    Some rows can include missing-data flags (e.g. -9999.0) under
    the VHELIO_AVG column. It means there is no enough information
    that star. To determine and filter, 'sential_abs_limit' is used.

    In addition, there may be no error values for APOGEE's VHELIO_AVG
    values. In such a case, the value of the apogee_manual_error parameter
    is used automatically as an error value. 
    """
    def _clean(col):
        arr = numpy.ma.filled(
            numpy.ma.masked_invalid(numpy.asarray(dataset[col], dtype=float)), fill_value=numpy.nan
            )
        arr = numpy.where(numpy.abs(arr) > sentinel_abs_limit, numpy.nan, arr)
        return arr

    gaia_radial_velocity = _clean(gaia_column)
    gaia_rv_err = _clean(gaia_error_column)
    apogee_radial_velocity = _clean(apogee_column)

    if apogee_error_column is not None:
        apogee_rv_err = _clean(apogee_error_column)
    else:
        apogee_rv_err = numpy.full(len(dataset), apogee_manual_error, dtype=float)

    missing_gaia = numpy.isnan(gaia_radial_velocity)
    combined_rv = numpy.where(missing_gaia, apogee_radial_velocity, gaia_radial_velocity)
    combined_rv_err = numpy.where(missing_gaia, apogee_rv_err, gaia_rv_err)
    still_missing = numpy.isnan(combined_rv)

    dataset = dataset.copy()
    dataset["radial_velocity_filled"] = Column(combined_rv, unit=units.km / units.s)
    dataset["radial_velocity_error_filled"] = Column(combined_rv_err, unit=units.km / units.s)
    dataset["radial_velocity_source"] = numpy.where(
        still_missing, "MISSING",
        numpy.where(missing_gaia, "APOGEE_VHELIO_AVG", "GAIA_DR3"),
    )

    n_total = len(dataset)
    n_gaia_missing = int(missing_gaia.sum())
    n_recovered = n_gaia_missing - int(still_missing.sum())
    print(
        f"Radial velocity: {n_gaia_missing}/{n_total} stars lacked a Gaia DR3 RV; "
        f"{n_recovered} recovered from APOGEE VHELIO_AVG; "
        f"{int(still_missing.sum())} still missing and will be dropped."
    )
    return dataset


#===============================================
# Cutting Radial Velocity Quality
#===============================================
##(See the Milky_Way_Mass_RotationCurve.ipynb notebook, Section 1.2.c)

def radial_velocity_quality_cut(dataset, rv_error_max, error_column="radial_velocity_error_filled"):
    error = dataset[error_column]
    mask = numpy.isfinite(error) & (error < rv_error_max)

    n_before = len(dataset)
    result = dataset[mask]
    print(f"radial_velocity_quality_cut: {n_before} -> {len(result)} stars (RV error < {rv_error_max} km/s)")

    return result

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

