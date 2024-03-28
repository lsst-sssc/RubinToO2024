import pytest

from numpy.testing import assert_allclose

from etc import get_exptime

class TestEtc:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.m5s = [23.70, 24.97, 24.52, 24.13, 23.56, 22.55]
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
