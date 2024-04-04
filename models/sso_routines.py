# author: Tim Lister <tlister@lco.global>
from datetime import datetime, timedelta

import requests
from astropy.time import Time
from astroquery.jplhorizons import Horizons

def get_horizons_ephem(rock, delta_days=15, ephem_step_size = '1h'):
    """
    Queries JPL HORIZONS for the ephemeris of a single object from the
    Simonyi Survey Telescope, Rubin Observatory [MPC code: X05] and produces 
    an ephemeris spanning close approach time +/- [delta_days] with steps
    of [ephem_step_size]. The default cutoffs are at airmass X=2.93 (~20d),
    no HA limit and no daylight cutoff.
  
    Parameters
    ----------
    rock list
        An entry from a query to the JPL SBDB close-approachers API. The only
        bits used are element [0] (object id) and element [2] (time of close
        approach, assumed to be a JD in the TDB timescale)

    Returns
    -------
    ephem `astropy.table.Table`
        ephemeris for the object. This is modified from the standard Horizons
        result returned by `astroquery` with the addition of a 'datetime' and
        'time_to_ca' (time to close approach) columns
    """

    # Hour angle limit (hours)
    ha_limit = 12
    # Airmass limit (corresponds to ~20d altitude)
    airmass_limit = 2.93
    # Obviously this should be True really but small number statistics of small
    # NEOs found in the S.H. means we don't care at this stage
    should_skip_daylight = False

    # Horizons quantities to return
    horizons_quantities = '1,3,4,9,19,20,23,24,38,42,33,46,47,48'

    # MPC Site Code (Rubin)
    site_code = 'X05'

    # Convert time of closest approach (str) to an Astropy Time object
    t_ca_tdb = Time(rock[2], format='jd', scale='tdb')
    # Convert to a UTC datetime and truncate to the minute
    t_ca_utc = t_ca_tdb.utc.datetime
    t_ca_utc_rounded = t_ca_utc.replace(second=0, microsecond=0)
    start_utc = t_ca_utc_rounded - timedelta(days=delta_days)
    end_utc   = t_ca_utc_rounded + timedelta(days=delta_days)

    # Initialize instance of class
    eph = Horizons(id=rock[0], id_type='smallbody', epochs={'start' : start_utc.strftime("%Y-%m-%d %H:%M"),
                'stop' : end_utc.strftime("%Y-%m-%d %H:%M"), 'step' : ephem_step_size}, location=site_code)

    ephem = eph.ephemerides(quantities=horizons_quantities,
        skip_daylight=should_skip_daylight, airmass_lessthan=airmass_limit,
        max_hour_angle=ha_limit)
    # Add datetime column
    dates = Time([datetime.strptime(d, "%Y-%b-%d %H:%M") for d in ephem['datetime_str']])
    if 'datetime' not in ephem.colnames:
        ephem.add_column(dates, name='datetime')
    # Make a time before/after close approach column. We just insert the time
    # difference as a float to prevent "ValueError: setting an array element with a sequence."
    # errors from the numpy maskedarray bug when e.g. plotting.
    ca_times = dates - t_ca_tdb.utc
    if 'time_to_ca' not in ephem.colnames:
        ephem.add_column(ca_times.value, name='time_to_ca')

    return ephem

def query_close_approachers(miss_distance='1LD', max_H=28.0, date_min=datetime(2023, 3, 14), date_max=datetime(2024, 3, 14)):
    """
    Queries the JPL Small Bodies Database Close-Approach Data API to find
    NEOs that meet the criteria passed as parameters
    API Docs: https://ssd-api.jpl.nasa.gov/doc/cad.html

    Parameters
    ----------
    miss_distance : `str`, optional
        The maximum value of the close-approach distance that will be returned.
        Default units for the JPL query are au; 'LD' can be appended for lunar
        distance. Default 1LD.
    max_H : `float` optional
        The maximum value of the absolute magnitude (H) that will be returned.
        Default 28 (~10 meters)
    date_min : `datetime`, optional
        The earliest time of close-approach that will be returned. Default 2023-03-14
    date_max : `datetime`, optional
        The latest time of close-approach that will be returned. Default 2024-03-14

    Returns
    -------
    close_approachers : (`dict`)
       A dictionary with four keys namely:
        "signature" : signature string containing API version,
        "count" : number of objects returned by the query (int),
        "fields": a list of fields for the data ["des","orbit_id","jd","cd","dist","dist_min","dist_max","v_rel","v_inf","t_sigma_f","h","diameter","diameter_sigma"],
        "data": a list of length <count> with data for each close approaching object,
            sorted by close approach distance
            Almost all fields, including numerical values, are encoded as strings
    """
    close_approachers = {}
    date_format = "%Y-%m-%d"
    query_url = f"https://ssd-api.jpl.nasa.gov/cad.api?dist-max={miss_distance}&date-min={date_min.strftime(date_format)}&date-max={date_max.strftime(date_format)}&h-max={max_H}&sort=dist&diameter=true"
    resp = requests.get(query_url)
    if resp.status_code == 200:
        close_approachers = resp.json()

    return close_approachers

def transform_Vmag(mag_V, passband, taxonomy='Mean'):
    """
    Returns the magnitude in <passband> for an asteroid with a V magnitude of
    <mag_V> and a taxonomic class of [taxonomy]. If the taxonomy is not given,
    a 'Mean' is assumed
    Taxonomy can be one of:
    'solar' - assuming solar colors (used by the MPC?),
    'mean'  - average of the S- and C-types is used,
    'neo'   - average weighted by the occurrence fraction among NEOs,
    's', 'c', 'q', 'x' - individual taxonomies

    Table 2. Asteroid magnitude transformations from Pan-STARRS1 AB filter magnitudes to the
    Johnson-Cousin V system based on Veres et al. (2015). Solar colors are also included for
    reference.
    Taxonomy    V-gP1   V-rP1   V-iP1   V-zP1   V-yP1   V-wP1
    Sun         -0.217  0.183   0.293   0.311   0.311   0.114
    Q           -0.312  0.252   0.379   0.238   0.158   0.156
    S           -0.325  0.275   0.470   0.416   0.411   0.199
    C           -0.238  0.194   0.308   0.320   0.316   0.120
    D           -0.281  0.246   0.460   0.551   0.627   0.191
    X           -0.247  0.207   0.367   0.419   0.450   0.146

    Mean (S+C)   -0.28  0.23    0.39    0.37    0.36    0.16

    According to Binzel et al. in _Asteroids IV_, p. 246:
    "About 90% of the known NEOs fall in the S-, Q-, C- and X-complexes
    with S- (50%), C- (15%), X- (10%) and Q- (10%) types dominating."
    """

    mag_mapping = { 'SOLAR' : {'g' : -0.217, 'r' : 0.183, 'i' : 0.293, 'z' : 0.311, 'Y' : 0.311, 'w' : 0.114},
                    'MEAN'  : {'g' : -0.28 , 'r' : 0.230, 'i' : 0.390, 'z' : 0.37 , 'Y' : 0.36 , 'w' : 0.160},
                       'S'  : {'g' : -0.325, 'r' : 0.275, 'i' : 0.470, 'z' : 0.416, 'Y' : 0.411, 'w' : 0.199},
                       'C'  : {'g' : -0.238, 'r' : 0.194, 'i' : 0.308, 'z' : 0.320, 'Y' : 0.316, 'w' : 0.120},
                       'Q'  : {'g' : -0.312, 'r' : 0.252, 'i' : 0.379, 'z' : 0.238, 'Y' : 0.158, 'w' : 0.156},
                       'X'  : {'g' : -0.247, 'r' : 0.207, 'i' : 0.367, 'z' : 0.419, 'Y' : 0.450, 'w' : 0.146},
                       'D'  : {'g' : -0.281, 'r' : 0.246, 'i' : 0.460, 'z' : 0.551, 'Y' : 0.627, 'w' : 0.191},
                     'NEO'  : {'g' : -0.254, 'r' : 0.213, 'i' : 0.356, 'z' : 0.322, 'Y' : 0.314, 'w' : 0.148},
                  }

    # Lookup taxonomy to try and get V-<passband> color terms
    color_terms = mag_mapping.get(taxonomy.upper(), None)

    # If we got a successful taxonomy lookup, try to lookup the <passband>
    # in the color terms
    delta_mag = None
    if color_terms:
        delta_mag = color_terms.get(passband, None)

    new_mag = None
    if delta_mag:
        new_mag = mag_V - delta_mag

    return new_mag
