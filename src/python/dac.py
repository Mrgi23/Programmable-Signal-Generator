import numpy as np
import scipy.signal as signal
from fir_filter import InverseSinc

class DAC():
    def __init__(self) -> None:
        self.inverse_sinc = InverseSinc()

    def __call__(
        self,
        digital: np.ndarray,
        mode: str,
        n_nyquist: int = 4,
        F_pass: float = 0.4,
        error_dB: float = 0.025
    ) -> np.ndarray:
        # Compute the reconstruction kernel.
        K = self.__kernel(mode, n_nyquist)

        # Filter the input signal with the sinc compensation filter, if necessary.
        if mode == "NRZ":
            b = self.inverse_sinc(F_pass, error_dB)
            filtered_digital = signal.lfilter(b, 1, digital)
        else:
            filtered_digital = digital

        # Define the analog signal.
        analog = np.zeros((len(digital) - 1) * n_nyquist + 1)

        # Upsample the input signal.
        analog[::n_nyquist] = filtered_digital[:]

        # Reconstruct the digital signal as analog.
        analog = np.convolve(K, analog)
        return analog

    def __kernel(self, mode: str, n_nyquist: int) -> np.ndarray:
        # Define kernel based on reconstruction function.
        if mode == "NRZ":
            # Zero-order hold.
            K = np.ones(n_nyquist)
            return K
        if mode == "RF":
            if n_nyquist % 2:
                raise ValueError("DAC: Indalid number of Nyquist zones for the RF mode.")

            # Bipolar zero-order hold.
            K = np.ones(n_nyquist)
            K[n_nyquist//2:] = -1
            return K
        raise ValueError("DAC: Invalid reconstruction mode.")
