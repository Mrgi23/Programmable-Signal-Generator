import numpy as np
import pytest
from interpolator import Interpolator

@pytest.fixture
def interpolator():
    return Interpolator()

def test_interpolator_valid_output(interpolator):
    f_max = 200
    fs = 1000
    t = np.linspace(0, 1, fs, endpoint=False)
    input = np.sin(2 * np.pi * f_max * t)

    output = interpolator(60.0, f_max, fs, input)
    output = np.abs(np.fft.fft(output))[:len(output)//2]
    assert(output.size == 16 * input.size // 2), "Invalid size of the interpolated signal."
    assert(np.isclose(output[f_max], fs // 2, atol=1e-5)), "Invalid frequency component."
    assert(np.all(np.isclose(output[output != output[f_max]], 0, atol=1e-5))), "Undefined frequency component."

def test_interpolator_filter_valid_output(interpolator):
    b = np.array([0.25, 0.25, 0.25, 0.25])
    input = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    output_expected = np.array([3.5, 4.5, 5.5, 6.5, 5.5, 4.5, 3.5, 2.5]).astype(complex)

    output = interpolator.filter(b, input)
    assert(output.size == input.size), "Invalid size of the filtered signal."
    assert(np.allclose(output, output_expected, atol=1e-5)), "Invalid filtered signal value."

def test_interpolator_filter_invalid_input(interpolator):
    with pytest.raises(ValueError, match="Filter must have at least one coefficient."):
        input = np.array([1, 2, 3, 4, 5])
        b = np.array([])
        interpolator.filter(b, input)

def test_interpolator_upsample_valid_output(interpolator):
    n = 3
    input = np.array([1, 2, 3, 4, 5])
    output_expected = np.array([1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0, 5, 0, 0]).astype(complex)

    output = interpolator.upsample(n, input)
    assert(output.size == n * input.size), "Invalid size of the upsampled signal."
    assert(np.allclose(output, output_expected, atol=1e-5)), "Invalid upsampled signal value."

def test_interpolator_upsample_invalid_input(interpolator):
    with pytest.raises(ValueError, match="Upsample factor must be greater than 0."):
        input = np.array([1, 2, 3, 4, 5])
        interpolator.upsample(0, input)
