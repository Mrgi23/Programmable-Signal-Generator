import numpy as np
import pytest
from complex_mixer import ComplexMixer

# Define the test object.
@pytest.fixture
def complex_mixer():
    return ComplexMixer()

def test_complex_mixer_valid_output(complex_mixer):
    # Define the input.
    f_shift = 150.00
    fs = 1000.00
    f = 300.00
    t = np.linspace(0, 1, int(fs), endpoint=False)
    I = np.cos(2 * np.pi * f * t)
    Q = np.zeros_like(t)

    # Compute the result.
    I_out = complex_mixer(f_shift, fs, I, Q)
    I_out = np.abs(np.fft.fft(I_out))

    # Test the result.
    spec = [int(f - f_shift), int(f + f_shift)]

    assert(I_out.size == I.size), "Invalid size of the shifted signal."
    assert(np.allclose(I_out[spec], fs / 4, atol=5e-2)), "Spectral components must be at f - f_shift and f + f_shift."
    assert(np.allclose(np.delete(I_out[:I_out.size//2], spec), 0, atol=5e-2)), "Noise must be close to 0."

def test_complex_mixer_invalid_input(complex_mixer):
    with pytest.raises(ValueError, match="ComplexMixer.__call__: Sampling frequency must be positive."):
        complex_mixer(200.0, -1000, np.array([1]), np.array([1]))
