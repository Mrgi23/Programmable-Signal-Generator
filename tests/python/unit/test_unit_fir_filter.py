import numpy as np
import pytest
from fir_filter import HalfBand, InverseSinc

# Define the test object.
@pytest.fixture
def inversesinc_filter():
    return InverseSinc()

def test_inversesinc_valid_output(inversesinc_filter):
    # Compute the result.
    b = inversesinc_filter(0.1, 0.1)

    # Test the result.
    assert(np.all(b[:len(b)//2] == b[-1:-len(b)//2:-1])), "Filter must be symetric."

def test_inversesinc_order_increase(inversesinc_filter):
    # Compute the result.
    b_lower = inversesinc_filter(0.1, 0.1)
    b_higher = inversesinc_filter(0.2, 0.1)

    # Test the result.
    assert(len(b_lower) < len(b_higher)), "Higher bandpass equals higher filter order."

    # Compute the result.
    b_lower = inversesinc_filter(0.1, 0.1)
    b_higher = inversesinc_filter(0.1, 0.01)

    # Test the result.
    assert(len(b_lower) < len(b_higher)), "Lower error equals higher filter order."

def test_inversesinc_too_high_order(inversesinc_filter):
    with pytest.raises(ValueError, match="InverseSinc.__call__: Too high order of the filter."):
        inversesinc_filter(0.49, 0.0025)

def test_inversesinc_invalid_input(inversesinc_filter):
    with pytest.raises(ValueError, match="InverseSinc.__call__: Passband must lie between 0.0 and 0.5."):
        inversesinc_filter(-0.1, 0.1)
    with pytest.raises(ValueError, match="InverseSinc.__call__: Passband must lie between 0.0 and 0.5."):
        inversesinc_filter(0.6, 0.1)

# Define the test object.
@pytest.fixture
def halfband_filter():
    return HalfBand()

def test_halfband_valid_output(halfband_filter):
    # Compute the result.
    coeffs = halfband_filter(120.0, 0.2)

    # Test the result.
    assert(np.isclose(coeffs[len(coeffs)//2], 0.5, atol=5e-3)), "Center coefficient should be close to 0.5."
    assert(np.all(coeffs[1::2][coeffs[1::2] != coeffs[len(coeffs)//2]]) == 0), "Every odd coefficient should be 0, except center coefficient."
    assert(np.all(coeffs[:len(coeffs)//2] == coeffs[-1:-len(coeffs)//2:-1])), "Filter must be symetric."

def test_halfband_order_increase(halfband_filter):
    # Compute the result.
    coeffs_lower = halfband_filter(120.0, 0.1)
    coeffs_higher = halfband_filter(120.0, 0.2)

    # Test the result.
    assert(len(coeffs_lower) < len(coeffs_higher)), "Higher bandpass equals higher filter order."

    # Compute the result.
    coeffs_lower = halfband_filter(120.0, 0.2)
    coeffs_higher = halfband_filter(130.0, 0.2)

    # Test the result.
    assert(len(coeffs_lower) < len(coeffs_higher)), "Higher attenuation equals higher filter order."

def test_halfband_too_high_order(halfband_filter):
    with pytest.raises(ValueError, match="HalfBand.__call__: Too high order of the filter."):
        halfband_filter(270.0, 0.24)

def test_halfband_invalid_input(halfband_filter):
    with pytest.raises(ValueError, match="Passband must lie between 0.0 and 0.25."):
        halfband_filter(120.0, -0.1)
    with pytest.raises(ValueError, match="Passband must lie between 0.0 and 0.25."):
        halfband_filter(120.0, 0.3)
