import numpy as np
import scipy.signal as signal
from abc import ABC

class FIR(ABC):
    """
    Base class for FIR-based filter designers.

    This abstract class stores the number of frequency-domain evaluation
    points used by derived FIR design algorithms. It does not implement
    filtering itself but provides shared configuration for subclasses.

    Parameters
    ----------
    n_points : int, optional
        Number of frequency-domain evaluation points used when validating
        filter responses.
    """
    def __init__(self, n_points: int = 8192) -> None:
        """
        Initialize the FIR.

        Parameters
        ----------
        n_points : int, optional
            Number of frequency-domain evaluation points used when validating
            filter responses.
        """
        super().__init__()
        self._n_points = n_points


class InverseSinc(FIR):
    """
    Inverse-sinc FIR filter designer.

    Designs an FIR filter that compensates for the sinc-shaped frequency
    response of zero-order hold (ZOH) DAC reconstruction. The filter is
    generated iteratively using least-squares FIR design until the passband
    error meets the specified tolerance.

    Parameters
    ----------
    n_points : int, optional
        Number of frequency-domain evaluation points used when validating
        the filter response.
    """

    def __init__(self, n_points: int = 8192) -> None:
        """
        Initialize the InverseSinc.

        Parameters
        ----------
        n_points : int, optional
            Number of evaluation points for frequency-response checks.
        """
        super().__init__(n_points)

    def __call__(self, F_pass: float, error_dB: float, n_spec: int = 16):
        """
        Design an inverse-sinc FIR filter.

        Generates a filter that approximates the inverse of the sinc response
        over the passband ``[0, F_pass]``. The design iteratively increases
        the filter order until the magnitude response satisfies the allowed
        ripple specified by ``error_dB``.

        Parameters
        ----------
        F_pass : float
            Passband edge. Must satisfy ``0 < F_pass < 0.5``.
        error_dB : float
            Allowed passband ripple in decibels.
        n_spec : int, optional
            Number of specification points used to construct the target
            frequency response.

        Returns
        -------
        ndarray
            FIR filter coefficients.

        Raises
        ------
        ValueError
            If ``F_pass`` is outside the valid range.
            If the required filter order exceeds internal limits.
        """
        if F_pass <= 0.0 or F_pass >= 0.5:
            raise ValueError("InverseSinc.__call__: Passband must lie between 0.0 and 0.5.")

        f = F_pass * np.linspace(0, 1, n_spec + 1)

        f_target = np.zeros(2 * (n_spec + 1))
        f_target[0::2] = 2 * f
        f_target[1::2] = 2 * (f + 1e-3)

        with np.errstate(divide='ignore', invalid='ignore'):
            h_target = np.zeros(2 * (n_spec + 1))
            h_target[0::2] = (np.pi * f) / np.sin(np.pi * f)
            h_target[0] = 1.0
            h_target[1::2] = (np.pi * f) / np.sin(np.pi * f)
            h_target[1] = 1.0

        N = 1

        while True:
            b = signal.firls(N, f_target, h_target)
            f_resp, h_inverse = signal.freqz(b, 1, worN=self._n_points, fs=1.0)

            with np.errstate(divide='ignore', invalid='ignore'):
                h_sinc = np.sin(f_resp * np.pi) / (f_resp * np.pi)
                h_sinc[0] = 1.0

            error = 10 ** (abs(error_dB) / 20)
            h_error = abs(h_sinc * h_inverse)[f_resp < F_pass]

            if np.all((h_error > 1 / error) & (h_error < error)):
                return b

            N += 2
            if N > 200:
                raise ValueError("InverseSinc.__call__: Too high order of the filter.")


class HalfBand(FIR):
    """
    Halfband low-pass FIR filter designer.

    Designs a halfband filter with a narrow transition band and symmetric
    coefficients. The filter is generated using the Parks-McClellan algorithm
    and validated against the desired passband ripple. If the design does not
    meet specifications, the filter order is increased iteratively.

    Parameters
    ----------
    n_points : int, optional
        Number of frequency-domain evaluation points used when validating
        the filter response.
    """

    def __init__(self, n_points: int = 8192) -> None:
        """
        Initialize the HalfBand.

        Parameters
        ----------
        n_points : int, optional
            Number of evaluation points for frequency-response checks.
        """
        super().__init__(n_points)

    def __call__(self, A_dB: float, F_pass: float) -> np.ndarray:
        """
        Design a halfband FIR filter.

        Creates a halfband low-pass filter with passband edge ``F_pass`` and
        attenuation ``A_dB``. The filter is designed using the Parks-McClellan
        algorithm and validated by checking the passband ripple. If the design
        fails, the filter order is increased until constraints are met.

        Parameters
        ----------
        A_dB : float
            Desired stopband attenuation in decibels.
        F_pass : float
            Passband edge. Must satisfy ``0 < F_pass < 0.25``.

        Returns
        -------
        ndarray
            Halfband FIR filter coefficients.

        Raises
        ------
        ValueError
            If ``F_pass`` is outside the valid range.
            If the required filter order exceeds internal limits.
        """
        if F_pass <= 0.0 or F_pass >= 0.25:
            raise ValueError("HalfBand.__call__: Passband must lie between 0.0 and 0.25.")

        delta_pass = 10 ** (-abs(A_dB) / 20)
        N = int(abs(A_dB) / (46 * (0.5 - 2 * F_pass)))

        if N % 2:
            N += 1

        while True:
            try:
                b = signal.remez(N, [0.0, 2 * F_pass], [1.0])
                f_resp, h = signal.freqz(b, 1, worN=self._n_points, fs=1.0)

                error = 2 * delta_pass
                h_error = abs(abs(h) - 1.0)[f_resp < 2 * F_pass]

                if np.all(h_error < error):
                    coeffs = np.zeros(2 * N - 1)
                    coeffs[::2] = b
                    coeffs[N - 1] = 1.0
                    return 0.5 * coeffs

            except ValueError:
                pass

            N += 2
            if N > 200:
                raise ValueError("HalfBand.__call__: Too high order of the filter.")
