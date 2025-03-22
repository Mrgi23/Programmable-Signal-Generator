import numpy as np
import pytest
from dac import DAC

def test_dac():
    # Define the test object.
    dac = DAC()

    # Define the input.
    fs = 1000.0
    f = 50.0
    t = np.linspace(0, 1, int(fs), endpoint=False)
    digital = np.cos(2 * np.pi * f * t)
    n_nyquist = 4
    F_pass = 0.4
    error_dB = 0.025

    # Compute the result.
    analog = dac(digital, "NRZ", n_nyquist, F_pass, error_dB)
    analog = np.abs(np.fft.fft(analog))

    # Test the result.
    spec = [int(f)]
    for i in range(1, n_nyquist // 2 + 1):
        spec.append(int(i * fs - f))
        spec.append(int(i * fs + f))
    spec.pop()
    max_value = np.max(analog[:analog.size//2])

    assert(analog.size == n_nyquist * digital.size), "Invalid size of the analog signal."
    assert(analog[int(f)] == max_value), "Original spectral component must be at f."
    assert(np.all(analog[spec] > 0.01 * max_value)), "Spectral replicas must be above 1% of maximum value."
    assert(np.all(np.delete(analog[:analog.size//2], spec) < 0.006 * max_value)), "Noise must be below 0.6% of maximum value."

    # Compute the result.
    analog = dac(digital, "RF", n_nyquist, F_pass, error_dB)
    analog = np.abs(np.fft.fft(analog))

    # Test the result.
    max_value = np.max(analog[:analog.size//2])

    assert(analog.size == n_nyquist * digital.size)
    assert(analog[int(fs - f)] == max_value), "Original spectral component must be at fs - f."
    assert(np.all(analog[spec] > 0.001 * max_value)), "Spectral replicas must be above 0.1% of maximum value."
    assert(np.allclose(np.delete(analog[:analog.size//2], spec), 0, atol=1e-5)), "Noise must be close to 0."
