import numpy as np
import pytest
from fir_filter import HalfBand

@pytest.fixture
def halfband_filter():
    return HalfBand()

def test_halfband_valid_output(halfband_filter):
    coeffs = halfband_filter(60.0, 0.2)
    assert(len(coeffs) > 0), "Filter coiefficients should not be empty array."
    assert(np.isclose(coeffs[len(coeffs)//2], 0.5, atol=1e-5)), "Center coefficient should be close to 0.5."
    assert(np.all(coeffs[1::2][coeffs[1::2] != coeffs[len(coeffs)//2]]) == 0), "Every odd coefficient should be 0, except center coefficient."

def test_halfband_invalid_argument(halfband_filter):
    f_passes = [-0.1, 0.3]
    for f_pass in f_passes:
        with pytest.raises(ValueError, match="Passband must lie between 0.0 and 0.25."):
            halfband_filter(60.0, f_pass)

def test_halfband_order_increase(halfband_filter):
    coeffs_lower = halfband_filter(60.0, 0.1)
    coeffs_higher = halfband_filter(60.0, 0.2)
    assert(len(coeffs_lower) < len(coeffs_higher)), "Higher bandpass equals higher filter order."

    coeffs_lower = halfband_filter(60.0, 0.2)
    coeffs_higher = halfband_filter(80.0, 0.2)
    assert(len(coeffs_lower) < len(coeffs_higher)), "Higher attenuation equals higher filter order."

def test_halfband_too_high_order(halfband_filter):
    with pytest.raises(ValueError, match="Too high order of the filter."):
        halfband_filter(100.0, 0.23)
