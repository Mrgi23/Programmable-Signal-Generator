import numpy as np
import scipy.signal as signal
from fir_filter import InverseSinc

class DAC:
    """
    Digital-to-analog reconstruction using selectable output modes.

    This class implements a simple DAC reconstruction model that converts a
    discrete-time digital signal into an analog-like waveform. It supports two
    reconstruction modes:

    - **"NRZ"**: Standard zero-order hold (ZOH) reconstruction with optional
      inverse-sinc compensation.
    - **"RF"**: Bipolar ZOH reconstruction used for RF DAC architectures.

    The class internally uses an `InverseSinc` filter to compensate for the
    inherent sinc-shaped frequency response of ZOH reconstruction when operating
    in NRZ mode.

    Parameters
    ----------
    n_points : int, optional
        Number of frequency-domain evaluation points used by the inverse-sinc
        filter designer.
    """

    def __init__(self, n_points: int = 8192) -> None:
        """
        Initialize the DAC.

        Parameters
        ----------
        n_points : int, optional
            Number of evaluation points for inverse-sinc filter design.
        """
        self.__inverse_sinc = InverseSinc(n_points)

    def __call__(
        self,
        digital: np.ndarray,
        mode: str,
        n_nyquist: int = 4,
        F_pass: float = 0.4,
        error_dB: float = 0.025
    ) -> np.ndarray:
        """
        Reconstruct an analog-like waveform from a digital input signal.

        The reconstruction process consists of:
        1. Generating a reconstruction kernel based on the selected mode.
        2. Optionally applying inverse-sinc compensation (NRZ mode only).
        3. Upsampling the signal by inserting zeros between samples.
        4. Convolving the upsampled signal with the reconstruction kernel.

        Parameters
        ----------
        digital : ndarray
            Input digital samples.
        mode : str
            Reconstruction mode. Must be `"NRZ"` or `"RF"`.
        n_nyquist : int, optional
            Number of Nyquist zones. Must be even when using `"RF"` mode.
        F_pass : float, optional
            Passband edge for inverse-sinc correction (NRZ mode).
        error_dB : float, optional
            Allowed ripple for inverse-sinc correction (NRZ mode).

        Returns
        -------
        ndarray
            Reconstructed analog-like output waveform.

        Raises
        ------
        ValueError
            If the reconstruction mode is invalid.
            If RF mode is selected with an odd number of Nyquist zones.
        """
        K = self.__kernel(mode, n_nyquist)

        if mode == "NRZ":
            b = self.__inverse_sinc(F_pass, error_dB)
            filtered_digital = signal.lfilter(b, 1, digital)
        else:
            filtered_digital = digital

        analog = np.zeros((len(digital) - 1) * n_nyquist + 1)
        analog[::n_nyquist] = filtered_digital

        analog = np.convolve(K, analog)
        return analog

    def __kernel(self, mode: str, n_nyquist: int) -> np.ndarray:
        """
        Generate the reconstruction kernel for the selected mode.

        The kernel defines how each digital sample is held or shaped during
        reconstruction:

        - **NRZ**: All-ones kernel (zero-order hold).
        - **RF**: First half +1, second half -1 (bipolar ZOH).

        Parameters
        ----------
        mode : str
            Reconstruction mode ("NRZ" or "RF").
        n_nyquist : int
            Number of Nyquist zones. Must be even for RF mode.

        Returns
        -------
        ndarray
            Reconstruction kernel.

        Raises
        ------
        ValueError
            If the mode is invalid.
            If RF mode is selected with an odd number of Nyquist zones.
        """
        if mode == "NRZ":
            return np.ones(n_nyquist)

        if mode == "RF":
            if n_nyquist % 2:
                raise ValueError("DAC.__kernel: Invalid number of Nyquist zones for the RF mode.")
            K = np.ones(n_nyquist)
            K[n_nyquist // 2:] = -1
            return K

        raise ValueError("DAC.__kernel: Invalid reconstruction mode.")
