import numpy as np
import pytest
from dac import DAC

@pytest.fixture
def converter():
    return DAC()

def test_dac_valid_output(mocker, converter):
    # Define the input.
    fs = 1000.0
    f = 50.0
    t = np.linspace(0, 1, int(fs), endpoint=False)
    digital = np.cos(2 * np.pi * f * t)
    n_nyquist = 4
    F_pass = 0.4
    error_dB = 0.025

    # Mock the InverseSinc filter.
    inverse_sinc_result = [np.array([1])]
    inverse_sinc_mock = mocker.Mock(side_effect=inverse_sinc_result)
    converter.inverse_sinc = inverse_sinc_mock

    # Compute the result.
    analog = converter(digital, "NRZ", n_nyquist, F_pass, error_dB)
    analog = np.abs(np.fft.fft(analog))

    # Test the result.
    spec = [int(f)]
    for i in range(1, 5):
        spec.append(int(i * fs - f))
        spec.append(int(i * fs + f))
    spec.pop()
    assert(analog.size == n_nyquist * fs), "Invalid size of the analog signal."
    assert(np.all(analog[spec] > 0.01 * np.max(analog))), "Spectral replicas must be above 1% of maximum value."
    assert(np.all(np.delete(analog, spec) < 0.006 * np.max(analog))), "Noise must be below 0.6% of maximum value."

    # Compute the result.
    analog = converter(digital, "RF", n_nyquist, F_pass, error_dB)
    analog = np.abs(np.fft.fft(analog))

    # Test the result.
    assert(analog.size == n_nyquist * fs)
    assert(np.all(analog[spec] > 0.001 * np.max(analog))), "Spectral replicas must be above 0.1% of maximum value."
    assert(np.allclose(np.delete(analog, spec), 0, atol=1e-5)), "Noise must be close to 0."

def test_dac_invalid_input(converter):
    with pytest.raises(ValueError, match="DAC.__kernel: Invalid reconstruction mode."):
        converter([1], "RFF")

    with pytest.raises(ValueError, match="DAC.__kernel: Indalid number of Nyquist zones for the RF mode."):
        converter([1], "RF", 5)
