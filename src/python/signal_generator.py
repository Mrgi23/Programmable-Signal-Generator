import numpy as np
from interpolator import Interpolator
from complex_mixer import ComplexMixer
from dac import DAC

class SignalGenerator:
    """
    Full digital signal-generation chain combining interpolation, frequency shifting,
    and DAC reconstruction.

    This class implements a complete multistage DSP pipeline that transforms a
    complex baseband signal into an analog-like waveform. The processing chain
    consists of:

    - **Interpolator**: Multistage x2 upsampling using halfband FIR filters.
    - **ComplexMixer**: Frequency translation using a CORDIC-based mixer.
    - **DAC**: Digital-to-analog reconstruction using NRZ or RF modes with optional
      inverse-sinc compensation.

    Parameters
    ----------
    n_steps : int, optional
        Number of x2 interpolation stages. The total interpolation factor is
        ``2**n_steps``.
    n_points : int, optional
        Number of frequency-domain evaluation points used by the filter designers.
    n_iter : int, optional
        Number of CORDIC iterations used by the complex mixer.
    f_res : float, optional
        Frequency resolution used by the complex mixer’s NCO.
    """

    def __init__(self, n_steps: int = 4, n_points: int = 8192, n_iter: int = 13, f_res: float = 1.0) -> None:
        """
        Initialize the signal-generation pipeline.

        Creates and configures the interpolator, complex mixer, and DAC components
        that form the full processing chain.

        Parameters
        ----------
        n_steps : int, optional
            Number of x2 interpolation stages.
        n_points : int, optional
            Number of evaluation points for filter design.
        n_iter : int, optional
            Number of CORDIC iterations for the complex mixer.
        f_res : float, optional
            Frequency resolution for the complex mixer’s NCO.
        """
        self.__interpolator = Interpolator(n_steps, n_points)
        self.__complex_mixer = ComplexMixer(n_iter, f_res)
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
        """
        Generate an analog-like waveform from a complex baseband signal.

        The processing chain is executed in three stages:

        1. **Interpolation**
           The input signal is upsampled by ``2**n_steps`` using halfband filters
           designed with attenuation ``A_dB`` and passband edge derived from
           ``f_max`` and ``fs``.

        2. **Frequency shifting**
           The interpolated I/Q components are mixed with a complex exponential
           at frequency ``f_shift`` using a CORDIC-based mixer.

        3. **DAC reconstruction**
           The shifted signal is converted into an analog-like waveform using the
           selected DAC mode (``"NRZ"`` or ``"RF"``), with optional inverse-sinc
           correction.

        Parameters
        ----------
        signal : ndarray
            Complex baseband input samples.
        fs : float
            Original sampling frequency. Must be strictly positive.
        f_max : float
            Maximum signal frequency present in the input.
        f_shift : float
            Frequency shift to apply after interpolation.
        mode : str
            DAC reconstruction mode (``"NRZ"`` or ``"RF"``).
        A_dB : float, optional
            Stopband attenuation for halfband filters.
        n_nyquist : int, optional
            Number of Nyquist zones for DAC reconstruction.
        F_pass : float, optional
            Passband edge for inverse-sinc correction.
        error_dB : float, optional
            Allowed ripple for inverse-sinc correction.

        Returns
        -------
        ndarray
            Final reconstructed analog-like waveform.
        """
        interpolated = self.__interpolator(A_dB, f_max, fs, signal)
        scale = 2 ** self.__interpolator.n_steps

        I = interpolated.real
        Q = interpolated.imag
        shifted = self.__complex_mixer(f_shift, scale * fs, I, Q)

        analog = self.__dac(shifted, mode, n_nyquist, F_pass, error_dB)
        return analog
