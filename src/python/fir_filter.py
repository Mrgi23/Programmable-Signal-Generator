from abc import ABC
import numpy as np
import scipy.signal as signal

class FIR(ABC):
    def __init__(self, n_points: int = 8192) -> None:
        super().__init__()
        self._n_points = n_points

class InverseSinc(FIR):
    def __init__(self, n_points = 8192) -> None:
        super().__init__(n_points)

    def __call__(self, F_pass: float, error_dB: float, n_spec: int = 16):
        if F_pass <= 0.0 or F_pass >= 0.5:
            raise ValueError("InverseSinc.__call__: Passband must lie between 0.0 and 0.5.")

        # Compute spectrum of frequencies.
        f = F_pass * np.linspace(0, 1, n_spec + 1)

        # Compute target frequencies.
        f_target = np.zeros(2 * (n_spec + 1))
        f_target[0::2] = 2 * f
        f_target[1::2] = 2 * (f + 1e-3)

        # Compute target frequency response (inverse sinc) and reduce singularity.
        with np.errstate(divide='ignore', invalid='ignore'):
            h_target = np.zeros(2 * (n_spec + 1))
            h_target[0::2] = (np.pi * f) / np.sin(np.pi * f)
            h_target[0] = 1.0
            h_target[1::2] = (np.pi * f) / np.sin(np.pi * f)
            h_target[1] = 1.0

        # Initial filter order. Type I filter, odd number of taps.
        N = 1

        # Iterate until filter is created.
        while True:
            # Design the filter.
            b = signal.firls(N, f_target, h_target)

            # Compute the frequency response.
            f, h_inverse = signal.freqz(b, 1, worN=self._n_points, fs=1.0)

            # Compute sinc response and reduce singularity.
            with np.errstate(divide='ignore', invalid='ignore'):
                h_sinc = np.sin(f * np.pi) / (f * np.pi)
                h_sinc[0] = 1.0

            # Calculate error and error frequency response.
            error = 10 ** (abs(error_dB) / 20)
            h_error = (abs(h_sinc * h_inverse)[f < F_pass])

            # Check if the filter satisfy constrains.
            if np.all((h_error > 1 / error) & (h_error < error)):
                return b

            # Increase the filter order.
            N += 2

            # Failsafe
            if N > 200:
                raise ValueError("InverseSinc.__call__: Too high order of the filter.")

class HalfBand(FIR):
    def __init__(self, n_points = 8192) -> None:
        super().__init__(n_points)

    def __call__(self, A_dB: float, F_pass: float) -> np.ndarray:
        if F_pass <= 0.0 or F_pass >= 0.25:
            raise ValueError("HalfBand.__call__: Passband must lie between 0.0 and 0.25.")

        # Passband ripple.
        delta_pass = 10 ** (-abs(A_dB) / 20)

        # Harris formula for initial filter order.
        N = int(abs(A_dB) / (46 * (0.5 - 2 * F_pass)))

        # Type II filter, even number of taps.
        if N % 2:
            N += 1

        # Iterate until filter is created.
        while True:
            try:
                # Design the filter.
                b = signal.remez(N, [0.0, 2 * F_pass], [1.0])

                # Compute the frequency response.
                f, h = signal.freqz(b, 1, worN=self._n_points, fs=1.0)

                # Calculate error and error frequency response.
                error = 2 * delta_pass
                h_error = abs(abs(h) - 1.0)[f < 2 * F_pass]

                # Check if the filter satisfy constrains.
                if np.all(h_error < error):
                    # Compute the halfband filter.
                    coeffs = np.zeros(2 * N - 1)
                    coeffs[::2] = b
                    coeffs[N-1] = 1.0
                    coeffs = 0.5 * coeffs
                    return coeffs
            except ValueError:
                pass

            # Increase the filter order.
            N += 2

            # Failsafe
            if N > 200:
                raise ValueError("HalfBand.__call__: Too high order of the filter.")
