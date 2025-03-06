import numpy as np
import pytest
from dac import DAC

@pytest.fixture
def converter():
    return DAC()

def test_dac_valid_output(converter):
    fs = 1000
    f = 50
    t = np.linspace(0, 1, fs, endpoint=False)
    digital = np.cos(2 * np.pi * f * t)
    n_nyquist = 4
    F_pass = 0.4
    error_dB = 0.025

    analog = converter(digital, "NRZ", n_nyquist, F_pass, error_dB)
    assert(analog.size == n_nyquist * fs), "Invalid size of the analog signal."
    analog = np.abs(np.fft.fft(analog))
    analog = analog[:(n_nyquist * fs)//2]
    assert(analog[f] == np.max(analog)), "Invalid analog signal value."
    assert(np.all(analog[[fs-f, fs+f, 2*fs-f]] > 0.01 * analog[f])), "Invalid analog signal value."
    assert(
        np.max(analog[
            (analog != analog[f]) &
            (analog != analog[fs-f]) &
            (analog != analog[fs+f]) &
            (analog != analog[2*fs-f])
        ]) < 0.006 * analog[f]
    ), "Invalid analog signal value."

    analog = converter(digital, "RF", n_nyquist, F_pass, error_dB)
    assert(analog.size == n_nyquist * fs)
    analog = np.abs(np.fft.fft(analog))
    analog = analog[:(n_nyquist * fs)//2]
    assert(analog[fs-f] == np.max(analog))
    assert(np.all(analog[[f, fs+f, 2*fs-f]] > 0.001 * analog[fs-f])), "Invalid analog signal value."
    assert(
        np.allclose(analog[
            (analog != analog[f]) &
            (analog != analog[fs-f]) &
            (analog != analog[fs+f]) &
            (analog != analog[2*fs-f])
        ], 0, atol=1e-5)
    ), "Invalid analog signal value."

def test_dac_invalid_input(converter):
    with pytest.raises(ValueError, match="DAC: Invalid reconstruction mode."):
        converter([1], "RFF")

    with pytest.raises(ValueError, match="DAC: Indalid number of Nyquist zones for the RF mode."):
        converter([1], "RF", 5)
