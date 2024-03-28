import pytest

import numpy as np
from numpy.testing import assert_allclose
import astropy.units as u

from etc import get_exptime, get_m5, calc_trailing_losses

class TestEtc:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.m5s = [23.70, 24.97, 24.52, 24.13, 23.56, 22.55]
        self.exptime = 30.0
        self.filters_all = "ugrizy"

    def test_darksky_exptime_zenith(self):
        expected_exptimes = { 'u': 30.0954,
                              'g': 30.0582,
                              'r': 30.2748,
                              'i': 30.2765,
                              'z': 30.1380,
                              'y': 29.9607
                            }

        for m5, filt in zip(self.m5s, self.filters_all):
            exptime_out = get_exptime(m5, filt, X=1.)
            assert_allclose(exptime_out, expected_exptimes[filt], rtol=1e-4)

    def test_darksky_exptime_am12(self):
        expected_exptimes = { 'u': 35.7850,
                              'g': 32.4761,
                              'r': 31.7601,
                              'i': 31.4127,
                              'z': 30.9253,
                              'y': 31.8972
                            }

        for m5, filt in zip(self.m5s, self.filters_all):
            exptime_out = get_exptime(m5, filt, X=1.2)
            assert_allclose(exptime_out, expected_exptimes[filt], rtol=1e-4)

    def test_twilight_exptime_zenith(self):
        expected_exptimes = { 'u': np.nan,
                              'g': np.nan,
                              'r': 442.1781,
                              'i': 541.1251,
                              'z': 473.6983,
                              'y': np.nan
                            }

        for m5, filt in zip(self.m5s, self.filters_all):
            exptime_out = get_exptime(m5, filt, X=1.0, twilight=True)
            assert_allclose(exptime_out, expected_exptimes[filt], rtol=1e-4)

    def test_twilight_exptime_alt20(self):
        expected_exptimes = { 'u': np.nan,
                              'g': np.nan,
                              'r': 700.2889,
                              'i': 770.7219,
                              'z': 606.7669,
                              'y': np.nan
                            }

        for m5, filt in zip(self.m5s, self.filters_all):
            exptime_out = get_exptime(m5, filt, X=2.92, twilight=True)
            assert_allclose(exptime_out, expected_exptimes[filt], rtol=1e-4)

    def test_darksky_m5_zenith(self):
        expected_m5s = { 'u': 23.6983,
                         'g': 24.9689,
                         'r': 24.5150,
                         'i': 24.1250,
                         'z': 23.5575,
                         'y': 22.5507
                       }

        for filt in self.filters_all:
            m5_out = get_m5(self.exptime, filt, X=1.)
            assert_allclose(m5_out, expected_m5s[filt], rtol=1e-4)

    def test_darksky_m5_am12(self):
        expected_m5s = { 'u': 23.6043,
                         'g': 24.9269,
                         'r': 24.4890,
                         'i': 24.1050,
                         'z': 23.5435,
                         'y': 22.5167
                       }

        for filt in self.filters_all:
            m5_out = get_m5(self.exptime, filt, X=1.2)
            assert_allclose(m5_out, expected_m5s[filt], rtol=1e-4)

    def test_twilight_m5_zenith(self):
        expected_m5s = { 'u': np.nan,
                         'g': np.nan,
                         'r': 23.0910,
                         'i': 22.5838,
                         'z': 22.0777,
                         'y': np.nan
                       }

        for filt in self.filters_all:
            m5_out = get_m5(self.exptime, filt, X=1.0, twilight=True)
            assert_allclose(m5_out, expected_m5s[filt], rtol=1e-4)

    def test_twilight_m5_am12(self):
        expected_m5s = { 'u': np.nan,
                         'g': np.nan,
                         'r': 23.0650,
                         'i': 22.5638,
                         'z': 22.0637,
                         'y': np.nan
                       }

        for filt in self.filters_all:
            m5_out = get_m5(self.exptime, filt, X=1.2, twilight=True)
            assert_allclose(m5_out, expected_m5s[filt], rtol=1e-4)

    def test_twilight_m5_alt20(self):
        expected_m5s = { 'u': np.nan,
                         'g': np.nan,
                         'r': 22.8414,
                         'i': 22.3918,
                         'z': 21.9433,
                         'y': np.nan
                       }

        for filt in self.filters_all:
            m5_out = get_m5(self.exptime, filt, X=2.92, twilight=True)
            assert_allclose(m5_out, expected_m5s[filt], rtol=1e-4)

class TestCalcTrailngLosses:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.exptime = 30.0
        self.max_rate = 10*u.deg/u.day

    def test_defaults(self):
        expected_losses = (0.5196914559265092, 0.8805600326623719)

        losses = calc_trailing_losses()

        assert_allclose(expected_losses, losses)

    def test_bad_floats(self):
        with pytest.raises(Exception):
            losses = calc_trailing_losses(2, 0.8, 15)

    def test_different_units(self):
        expected_losses = (0.5196914559265092, 0.8805600326623719)

        losses = calc_trailing_losses(5*u.arcsec/u.min, exptime=0.5*u.min)

        assert_allclose(expected_losses, losses)
