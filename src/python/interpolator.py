import numpy as np
import scipy.signal as signal
from fir_filter import HalfBand

class Interpolator():
    def __init__(self) -> None:
        self.halfband = HalfBand()

    def __call__(self, A_dB: float, f_max: float, fs: float, input: np.ndarray) -> np.ndarray:
        if fs <= 0:
            raise ValueError("Interpolator.__call__: Sampling frequency must be positive.")

        # Dedfine the output signal.
        output = input

        # Propagate output signal through the interpolation 4 times.
        for i in range(4):
            # Upsample previous output signal by factor 2.
            output = self.__upsample(2, output)

            # Compute FIR filter coefficients.
            factor = (2 ** i)
            F_pass = f_max / (factor * fs)
            b = self.halfband(A_dB, F_pass)

            # Filter the upsampled signal.
            output = self.__filter(b, output)
        return output

    def __filter(self, b: np.ndarray, input: np.ndarray) -> np.ndarray:
        # Prepare signal for filtering by expanding its size for len(b) elements.
        output = np.concatenate((input, input[:len(b)]))

        # Filter signal.
        output = signal.lfilter(b, 1, output)

        # Remove additional elements of the signal.
        output = output[len(b):]
        return output

    def __upsample(self, n: int, input: np.ndarray) -> np.ndarray:
        # Define empty output with the valid size.
        output = np.zeros(len(input) * n, dtype=complex)

        # Upsample signal by adding n-1 zeros between every element.
        output[::n] = input[:]
        return output
