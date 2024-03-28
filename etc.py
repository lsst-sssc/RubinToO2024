# author: Igor Andreoni <igor.andreoni@gmail.com.>
# author: Tim Lister <tlister@lco.global>
# Formulae from https://smtn-002.lsst.io/v/OPSIM-1134/index.html
# with delta magnitudes for the twilight survey and trailing loss code from Lynne Jones

import warnings
import numpy as np
from astropy import units as u

params = {
          "u": {"Cm": 22.97,
                "dCm_inf": 0.54,
                "zp": 26.52,
                "fwhm": 0.92,
                "m_darksky": 23.05,
                "k_atm": 0.47
                },
          "g": {"Cm": 24.58,
                "dCm_inf": 0.09,
                "zp": 28.51,
                "fwhm": 0.87,
                "m_darksky": 22.25,
                "k_atm": 0.21
                },
          "r": {"Cm": 24.6,
                "dCm_inf": 0.04,
                "zp": 28.36,
                "fwhm": 0.83,
                "fwhm_twilight": 1.43,
                "m_darksky": 21.2,
                "m_twilight": 19.47, # median value
                "k_atm": 0.13
                },
          "i": {"Cm": 24.54,
                "dCm_inf": 0.03,
                "zp": 28.17,
                "fwhm": 0.80,
                "fwhm_twilight": 1.49,
                "m_darksky": 20.46,
                "m_twilight": 18.68, # median value
                "k_atm": 0.10
                },
          "z": {"Cm": 24.37,
                "dCm_inf": 0.02,
                "zp": 27.78,
                "fwhm": 0.78,
                "fwhm_twilight": 1.42,
                "m_darksky": 19.61,
                "m_twilight": 17.92, # median value
                "k_atm": 0.07
                },
          "y": {"Cm": 23.84,
                "dCm_inf": 0.02,
                "zp": 26.82,
                "fwhm": 0.76,
                "m_darksky": 18.6,
                "k_atm": 0.17
                }
          }

@u.quantity_input(velocity=u.deg/u.day, seeing=u.arcsec, exptime=u.s)
def calc_trailing_losses(velocity=2*u.deg/u.day, seeing=0.8*u.arcsec, exptime=30.0*u.s):
    """Calculate the detection and SNR trailing losses.
    Code ported from the original Lynne Jones rubin_sim.moving_objects.base_obs
    code and modified slightly to require AstroPy units function arguments
    by Tim Lister.

    'Trailing' losses = loss in sensitivity due to the photons from the
    source being spread over more pixels; thus more sky background is
    included when calculating the flux from the object and thus the SNR
    is lower than for an equivalent brightness stationary/PSF-like source.
    dmag_trail represents this loss.

    'Detection' trailing losses = loss in sensitivity due to the photons
    from the source being spread over more pixels, in a non-stellar-PSF
    way, while source detection is (typically) done using a stellar PSF
    filter and 5-sigma cutoff values based on assuming peaks from
    stellar PSF's above the background; thus the SNR is lower than for an
    equivalent brightness stationary/PSF-like source (and by a greater
    factor than just the simple SNR trailing loss above).
    dmag_detect represents this loss.

    Parameters
    ----------
    velocity : `astropy.units.Quantity`
        The velocity of the moving objects, in deg/day.
    seeing : `astropy.units.Quantity`
        The seeing of the images, in arcseconds.
    exptime : `astropy.units.Quantity`, optional
        The exposure time of the images, in seconds. Default 30.

    Returns
    -------
    dmag trail, dmag_detect : (`np.ndarray` `np.ndarray`)
    or (`float`, `float`)
        dmag_trail and dmag_detect for each set of
        velocity/seeing/texp values.
    """

    a_trail = 0.761
    b_trail = 1.162
    a_det = 0.420
    b_det = 0.003
    x = velocity.to(u.deg/u.day).value * exptime.to(u.s).value / seeing.to(u.arcsec).value / 24.0
    dmag_trail = 1.25 * np.log10(1 + a_trail * x**2 / (1 + b_trail * x))
    dmag_detect = 1.25 * np.log10(1 + a_det * x**2 / (1 + b_det * x))
    return (dmag_trail, dmag_detect)

def calc_event_time_budget(n_fields=1, filters=['griz'], exptimes=[30, 30, 30, 30]):
    """Calculates the time to follow a single ToO event in the specified filters
    and exposure times. Assumes only a single repeat of the observation pattern
    and not a (repeated) cadence.

    Parameters
    ----------
    n_fields : `int`
        The number of fields to be tiled for the event
    filters : `list`
        The filters to be observed. If the full set of ['ugrizy'] are specified,
        one will be dropped and the exposure times for u and z will be scaled
        by the ratio of the mounted time.
    exptimes : `list`, optional
        The exposure time of the images, in each filter in seconds. Default 30.

    Returns
    -------
    budget_strategy_nfields_hr : `float`
        total time (in hours) to execute the strategy over the <n_fields>.
    """

    # Overheads
    overhead_between_exposures = 7.0  # in seconds

    # Overhead first slew per visit, assuming it may be further away than tiling
    overhead_first_slew = 30

    # Overhead for slew and settle, assuming closer fields beyond the initial slew (not currently used)
    overhead_slew_settle = 10

    # Overhead for filter change
    overhead_filter_change = 120.0  # in seconds

     # Initialize total overhead slews
    overhead_first_slew_total = 0
    # Initialize total overhead readouts
    overhead_filter_change_total = 0
    # Initialize total overhead readouts
    overhead_between_exposures_total = 0
    # Initialize total exposure time
    exptime_total = 0

    # Main calc
    # Add overheads for the first slew
    overhead_first_slew_total += overhead_first_slew
    # Initialize
    exptime_visit = 0
    # Check if 6 filters are given
    if len(filters) == 6:
        # add overhead between exposures
        overhead_between_exposures_visit = overhead_between_exposures * 4
        # Add overhead change of filter, assuming the first filter changes, too
        overhead_filter_change_visit =  overhead_filter_change * 5
        # If 6 filters are given, one must be dropped.
        # Indicate here the number of days per month in which u band is used.
        # Y band will be assumed for the other visits. The total time budget will be a combination of the two.
        u_days_ratio = 14./30
        # Combine exposure times between u and y
        for filt, exptime in zip(filters, exptimes):
            if filt == "u":
                exptime_visit += exptime * u_days_ratio
            elif filt == "y":
                exptime_visit += exptime * (1-u_days_ratio)
            else:
                exptime_visit += exptime
    # In case <= 5 filters are given for this visit
    else:
        overhead_between_exposures_visit = overhead_between_exposures * len(exptimes) # Not sure why the original '-1' - still have to pay final readout cost shirley?
        # Add overhead change of filter, assuming the first filter changes, too
        overhead_filter_change_visit = overhead_filter_change * len(exptimes)
        # Exposure time
        exptime_visit = np.sum(exptimes)
        # Add the exposure times to the epoch to the total
        exptime_total += exptime_visit
        # Add overheads
        overhead_between_exposures_total += overhead_between_exposures_visit
        overhead_filter_change_total += overhead_filter_change_visit
        # Print results for the epoch
        print(f"  Epoch T+4 hr:")
        print(f"    exposure times: {'{:.0f}'.format(exptime_visit)}s")
        print(f"    overhead change filter: {overhead_filter_change_visit}s")
        print(f"    overhead between exposures: {overhead_between_exposures_visit}s")

    # calculate the total overheads, convert overhead from seconds to hours
    # IMPORTANT: the filter change and first slew overhead must to be divided by the number of fields,
    # if we complete the tiling before changing filters
    overheads_total_nfields = np.sum((overhead_between_exposures_total * n_fields)/60/60 +
                       overhead_first_slew_total/n_fields/60/60 +
                       overhead_filter_change_total/n_fields/60/60)

    # calculate the total exposure time per field/pointing in hours
    total_exposure_time_hr = exptime_total/60/60

    # Calculate the total time budget in hours per event
    budget_strategy_nfields_hr = total_exposure_time_hr * n_fields + overheads_total_nfields

    print(f"Total exposure time per pointing: {'{:.3f}'.format(total_exposure_time_hr)}hr")
    #print(f"Total overheads per pointing, assuming {n_fields} fields: {'{:.2f}'.format(overheads_total_nfields)}hr ({'{:.2f}'.format(100*overheads/budget_strategy_hr)}% of total budget)")
    print(f"Total time for {n_fields} fields: {'{:.3f}'.format(budget_strategy_nfields_hr)}hr")
    print("-- \n")

    return budget_strategy_nfields_hr

def get_exptime(m5, filt, X=1.0, twilight=False):
    """
    Given a certain depth, return the exposure time

    Parameters
    ----------
    m5 float
        depth 5sigma (mag)
    filt str
        filter (one of ugrizy)
    X float
        airmass
    twilight bool
        Whether to use the twilight survey numbers

    Returns
    -------
    exptime float
        exposure time to reach limiting magnitude
    """
    # assign the parameters to variables
    Cm = params[filt]["Cm"]
    k_atm = params[filt]["k_atm"]
    fwhm = params[filt]["fwhm"]
    m_darksky = params[filt]["m_darksky"]
    # Important: assuming darksky
    m_sky = m_darksky
    if twilight:
        m_sky = params[filt].get('m_twilight', -99.0)
        fwhm = params[filt].get("fwhm_twilight", -99.0)
        # Suppress warnings for the invalid filters
        warnings.simplefilter('ignore', RuntimeWarning)
    # FIXME approximation dCm=0 (fine within 0.3s for 30s exposures)
    # dCm_inf = params[filt]["dCm_inf"]
    # Tscale = exptime / 30. * 10**(-1 * 0.4 * (m_sky - m_darksky))
    # dCm = dCm_inf - 1.25 * np.log10(1 + (10**(0.8 * dCm_inf) - 1) / Tscale)
    dCm = 0
    # Calculate the exposure time
    exptime = 30 * 10 ** ((1 / 1.25) * (m5 - Cm - dCm - 0.5 * (m_sky - 21.) -
                                        2.5 * np.log10(0.7 / fwhm) +
                                        k_atm*(X - 1.0)))

    return exptime


def get_m5(exptime, filt, X=1.0, twilight=False):
    """
    Given a certain exposure time return 5sigma depth
    Parameters
    ----------
    exptime int or float
        exposure time
    filt str
        filter (one of ugrizy)
    X float
        airmass
    twilight bool
        Whether to use the twilight survey numbers

    Returns
    -------
    m5 float
        5sigma limiting magnitude
    """
    # assign the parameters to variables
    Cm = params[filt]["Cm"]
    dCm_inf = params[filt]["dCm_inf"]
    k_atm = params[filt]["k_atm"]
    fwhm = params[filt]["fwhm"]
    m_darksky = params[filt]["m_darksky"]
    # Important: assuming darksky
    m_sky = m_darksky
    if twilight:
        m_sky = params[filt].get('m_twilight', -99.0)
        fwhm = params[filt].get("fwhm_twilight", -99.0)
        # Suppress warnings for the invalid filters
        warnings.simplefilter('ignore', RuntimeWarning)
    # Calculate m5
    Tscale = exptime / 30. * 10**(-1 * 0.4 * (m_sky - m_darksky))
    dCm = dCm_inf - 1.25 * np.log10(1 + (10**(0.8 * dCm_inf) - 1) / Tscale)
    m5 = Cm + dCm + 0.5 * (m_sky - 21.) + 2.5 * np.log10(0.7 / fwhm) + \
        1.25 * np.log10(exptime / 30.) - k_atm*(X - 1.0)

    return m5


if __name__ == "__main__":
    # Testing getting exptime given a list of mag limits
    m5s = [23.70, 24.97, 24.52, 24.13, 23.56, 22.55]
    for m5, filt in zip(m5s, "ugrizy"):
        # m5_out = get_m5(exptime, filt, X=1.)
        exptime_out = get_exptime(m5, filt, X=1.)
        # print('{:.2f}'.format(m5_out))
        print('{:.2f}'.format(exptime_out))
