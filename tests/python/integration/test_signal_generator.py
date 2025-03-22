import numpy as np
import pytest
from signal_generator import SignalGenerator

def test_signal_generator():
    # Define the test object.
    signal_generator = SignalGenerator()
    M = 2 ** signal_generator._SignalGenerator__interpolator.n_steps

    # Define the input.
    fs = 60
    f_min = 12
    f_max = 24
    f_shift_nrz = 3 / 8 * fs * M
    f_shift_rf = 6 / 8 * fs * M
    n_nyquist = 4

    data = np.loadtxt("./data/testSignal.txt", delimiter=",")
    signal = data[:, 0] + 1j * data[:, 1]
    N = signal.size

    # Compute the result.
    analog = signal_generator(signal, fs, f_max, f_shift_nrz, "NRZ", 60.0, n_nyquist)
    analog = np.abs(np.fft.fft(analog))

    # Test the result.
    spec_max = []
    spec = []
    f0 = f_shift_nrz + f_min
    f1 = f_shift_nrz + f_max
    for i in range(analog.size // 2):
        if (i >= int(f0 / fs * N) and i < int(f1 / fs * N)):
            spec_max.append(i)
        elif (
            (i > int((M * fs - f1) / fs * N) and i <= int((M * fs - f0) / fs * N)) or
            (i >= int((M * fs + f0) / fs * N) and i < int((M * fs + f1) / fs * N)) or
            (i > int((2 * M * fs - f1) / fs * N) and i <= int((2 * M * fs - f0) / fs * N))
        ):
            spec.append(i)
    max_value = np.max(analog[:analog.size//2])

    assert(analog.size == n_nyquist * M * N), "Invalid size of the generated signal."
    assert(np.allclose(analog[spec_max], max_value, atol=5e-2)), "Original spectral component must be between f_shift + f_min and f_shift + f_max."
    assert(np.all(analog[spec] > 0.01 * max_value)), "Spectral replicas must be above 1% of maximum value."
    assert(np.all(np.delete(analog[:analog.size//2], spec_max + spec) < 0.006 * max_value)), "Noise must be below 0.6% of maximum value."

    # Compute the result.
    analog = signal_generator(signal, fs, f_max, f_shift_rf, "RF", 60.0, n_nyquist)
    analog = np.abs(np.fft.fft(analog))

    # Test the result.
    spec_max = []
    spec = []
    f0 = f_shift_rf + f_min
    f1 = f_shift_rf + f_max
    for i in range(analog.size // 2):
        if (i >= int(f0 / fs * N) and i < int(f1 / fs * N)):
            spec_max.append(i)
        elif (
            (i > int((M * fs - f1) / fs * N) and i <= int((M * fs - f0) / fs * N)) or
            (i >= int((M * fs + f0) / fs * N) and i < int((M * fs + f1) / fs * N)) or
            (i > int((2 * M * fs - f1) / fs * N) and i <= int((2 * M * fs - f0) / fs * N))
        ):
            spec.append(i)
    max_value = np.max(analog[:analog.size//2])

    assert(analog.size == n_nyquist * M * N), "Invalid size of the generated signal."
    assert(np.allclose(analog[spec_max], max_value, atol=5e-2)), "Original spectral component must be between f_shift + f_min and f_shift + f_max."
    assert(np.all(analog[spec] > 0.001 * max_value)), "Spectral replicas must be above 0.1% of maximum value."
    assert(np.allclose(np.delete(analog[:analog.size//2], spec_max + spec), 0, atol=1e-3)), "Noise must be close to 0."
