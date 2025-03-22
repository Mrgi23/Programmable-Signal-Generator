import numpy as np
import pytest
from interpolator import Interpolator

def test_interpolator():
    # Define the test object.
    interpolator = Interpolator()

    # Define the input.
    f_max = 400.00
    fs = 1000.00
    t = np.linspace(0, 1, int(fs), endpoint=False)
    input = np.cos(2 * np.pi * f_max * t)

    # Compute the result.
    output = interpolator(120.0, f_max, fs, input)
    output = np.abs(np.fft.fft(output))
    N = 2 ** interpolator.n_steps

    # Test the result.
    assert(output.size == N * input.size), "Invalid size of the interpolated signal."
    assert(np.isclose(output[int(f_max)], fs / 2, atol=5e-3)), "Spectral component must be at fmax."
    assert(np.allclose(np.delete(output[:output.size//2], int(f_max)), 0, atol=5e-3)), "Noise must be close to 0."
