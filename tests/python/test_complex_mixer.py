import numpy as np
import pytest
from complex_mixer import ComplexMixer

@pytest.fixture
def complex_mixer():
    return ComplexMixer()

def test_complex_mixer_valid_output(complex_mixer):
    f_shift = 250
    fs = 1000
    f = 100
    t = np.linspace(0, 1, fs, endpoint=False)
    I = np.cos(2 * np.pi * f * t)
    Q = np.zeros_like(t)

    I_out = complex_mixer(f_shift, fs, I, Q)
    assert(I_out.size == I.size), "Invalid size of the shifted signal."
    I_out = np.abs(np.fft.fft(I_out))[:len(I_out)//2]
    assert(np.allclose(I_out[[f_shift-f, f_shift+f]], fs // 4, atol=1e-5)), "Spectral components must be at fshift-f and fshifted+f."
    assert(np.allclose(I_out[(I_out != I_out[f_shift-f]) & (I_out != I_out[f_shift+f])], 0, atol=1e-5)), "Noise must be close to 0."

def test_complex_mixer_invalid_input(complex_mixer):
    with pytest.raises(ValueError, match="ComplexMixer.__call__: Sampling frequency must be positive."):
        complex_mixer(200.0, -1000, np.array([1]), np.array([1]))
