import numpy as np
import pytest
from interpolator import Interpolator

@pytest.fixture
def interpolator():
    return Interpolator()

def test_interpolator_valid_output(mocker, interpolator):
    # Define the input.
    f_max = 400.00
    fs = 1000.00
    t = np.linspace(0, 1, int(fs), endpoint=False)
    input = np.sin(2 * np.pi * f_max * t)

    # Mock the HalfBand filter.
    halfband_result = [np.array([1])] * interpolator.N
    halfband_mock = mocker.Mock(side_effect=halfband_result)
    interpolator._Interpolator__halfband = halfband_mock

    # Compute the result.
    output = interpolator(120.0, f_max, fs, input)
    output = np.abs(np.fft.fft(output))

    # Test the result.
    spec = []
    for i in range(2 ** interpolator.N):
        spec.append(int(i * fs + f_max))
        spec.append(int((i + 1) * fs - f_max))
    assert(output.size == (2 ** interpolator.N) * input.size), "Invalid size of the interpolated signal."
    assert(np.allclose(output[spec], fs / 2, atol=1e-9)), "Spectral component must be at fmax."
    assert(np.allclose(np.delete(output, spec), 0, atol=1e-9)), "Noise must be close to 0."

def test_interpolator__invalid_input(interpolator):
    with pytest.raises(ValueError, match="Interpolator.__call__: Sampling frequency must be positive."):
        interpolator(60.0, 200.0, -1000.0, np.array([1]))
