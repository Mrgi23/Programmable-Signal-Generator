import numpy as np
import scipy.signal as signal
from fir_filter import HalfBand

class Interpolator:
    """
    Multistage complex-signal interpolator using halfband FIR filters.

    This class performs repeated x2 interpolation of a complex input signal.
    Each interpolation stage upsamples the signal by inserting zeros and then
    applies a halfband low-pass filter to suppress spectral images. The number
    of interpolation stages is configurable, and the halfband filter is designed
    adaptively based on the desired attenuation and passband edge.

    Parameters
    ----------
    n_steps : int, optional
        Number of x2 interpolation stages. The total interpolation factor is
        ``2**n_steps``.
    n_points : int, optional
        Number of frequency-domain evaluation points used by the halfband
        filter designer.
    """

    def __init__(self, n_steps: int = 4, n_points: int = 8192) -> None:
        """
        Initialize the Interpolator.

        Parameters
        ----------
        n_steps : int, optional
            Number of x2 interpolation stages.
        n_points : int, optional
            Number of evaluation points for halfband filter design.
        """
        self.__n_steps = n_steps
        self.__halfband = HalfBand(n_points)

    @property
    def n_steps(self):
        """
        Number of x2 interpolation stages.

        Returns
        -------
        int
            Number of interpolation stages.
        """
        return self.__n_steps

    def __call__(self, A_dB: float, f_max: float, fs: float, input: np.ndarray) -> np.ndarray:
        """
        Perform multistage interpolation on a complex input signal.

        The signal is repeatedly upsampled by a factor of two and filtered using
        a halfband low-pass filter. The passband edge is adjusted at each stage
        based on the maximum signal frequency ``f_max`` and the sampling
        frequency ``fs``. The halfband filter is designed using the specified
        attenuation ``A_dB``.

        Parameters
        ----------
        A_dB : float
            Desired stopband attenuation for the halfband filter.
        f_max : float
            Maximum signal frequency present in the input.
        fs : float
            Sampling frequency. Must be strictly positive.
        input : ndarray
            Complex input samples.

        Returns
        -------
        ndarray
            Interpolated complex output signal.

        Raises
        ------
        ValueError
            If ``fs <= 0``.
        """
        if fs <= 0:
            raise ValueError("Interpolator.__call__: Sampling frequency must be positive.")

        output = input

        for i in range(1, self.n_steps + 1):
            output = self.__upsample(2, output)
            factor = 2 ** i
            F_pass = f_max / (factor * fs)
            b = self.__halfband(A_dB, F_pass)
            output = self.__filter(b, output)

        return output

    def __filter(self, b: np.ndarray, input: np.ndarray) -> np.ndarray:
        """
        Apply an FIR filter to a complex signal.

        The input is padded to avoid circular convolution artifacts, filtered
        using direct-form FIR filtering, and then trimmed back to the correct
        length.

        Parameters
        ----------
        b : ndarray
            FIR filter coefficients.
        input : ndarray
            Complex input samples.

        Returns
        -------
        ndarray
            Filtered complex output signal.
        """
        output = np.concatenate((input, input[:len(b)]))
        output = signal.lfilter(b, 1, output)
        output = output[len(b):]
        return output

    def __upsample(self, n: int, input: np.ndarray) -> np.ndarray:
        """
        Upsample a complex signal by inserting zeros.

        Parameters
        ----------
        n : int
            Upsampling factor.
        input : ndarray
            Complex input samples.

        Returns
        -------
        ndarray
            Upsampled complex signal.
        """
        output = np.zeros(len(input) * n, dtype=complex)
        output[::n] = input
        return output
