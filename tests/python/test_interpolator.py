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
    assert(output.size == 16 * input.size), "Invalid size of the interpolated signal."
    output = np.abs(np.fft.fft(output))[:len(output)//2]
    assert(np.isclose(output[f_max], fs // 2, atol=1e-5)), "Spectral component must be at fmax."
    assert(np.allclose(output[output != output[f_max]], 0, atol=1e-5)), "Noise must be close to 0."

def test_interpolator__invalid_input(interpolator):
    with pytest.raises(ValueError, match="Interpolator.__call__: Sampling frequency must be positive."):
        interpolator(60.0, 200.0, -1000.0, np.array([1]))
