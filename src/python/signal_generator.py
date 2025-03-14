import numpy as np
from interpolator import Interpolator
from complex_mixer import ComplexMixer
from dac import DAC

class SignalGenerator():
    def __init__(self, N: int = 4, n_points: int = 8192, n_iter: int = 13, f_res: float = 1.0) -> None:
        # Define the interpolator.
        self.__interpolator = Interpolator(N, n_points)

        # Define the complex mixer.
        self.__complex_mixer = ComplexMixer(n_iter, f_res)

        # Define the digital-to-analog converter.
        self.__dac = DAC(n_points)

    def __call__(
        self,
        signal: np.ndarray,
        fs: float,
        f_max: float,
        f_shift: float,
        mode: str,
        A_dB: float = 60.0,
        n_nyquist: int = 4,
        F_pass: float = 0.4,
        error_dB: float = 0.025
    ) -> np.ndarray:
        # Interpolate the input signal.
        interpolated = self.__interpolator(A_dB, f_max, fs, signal)
        scale = len(interpolated) // len(signal)

        # Shift the interpolated signal.
        I = interpolated.real
        Q = interpolated.imag
        shifted = self.__complex_mixer(f_shift, scale * fs, I, Q)

        # Convert the digital signal to analog.
        analog = self.__dac(shifted, mode, n_nyquist, F_pass, error_dB)
        return analog
