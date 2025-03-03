from abc import ABC
import numpy as np
import scipy.signal as signal

class FIR(ABC):
    def __init__(self, n_points: int = 8192) -> None:
        super().__init__()
        self.n_points = n_points


class HalfBand(FIR):
    def __init__(self, n_points = 8192):
        super().__init__(n_points)

    def __call__(self, A_dB: float, F_pass: float) -> np.ndarray:
        # Specifications are not valid.
        if F_pass < 0.0 or F_pass > 0.25:
            raise ValueError("HalfBand.__call__: Passband must lie between 0.0 and 0.25.")

        # Passband ripple.
        delta_pass = 10 ** (-abs(A_dB) / 20)

        # Harris formula for initial filter order.
        N = int(2 * abs(A_dB) / (23 * (0.5 - 2 * F_pass)))

        # Force even number of taps, since scipy.remez uses numtaps instead of order.
        if N % 2:
            N += 1

        # Iterate until filter is created.
        while True:
            # Design the filter.
            try:
                b = signal.remez(N, [0.0, 2*F_pass], [1.0], weight=[1.0])
            except ValueError:
                b = np.zeros(N)
            # Calculate the frequency response & the amplitude
            w, h = signal.freqz(b, 1, worN=self.n_points, fs=1)
            H = abs(h)

            if not (np.sum((w < 2 * F_pass) * abs(H - 1.0) > 2 * delta_pass)):
                # Create the halfband filter
                coeffs = np.zeros(2*N-1)
                coeffs[::2] = b
                coeffs[N-1] = 1.0
                coeffs = 0.5 * coeffs
                return coeffs

            # Increase the filter order.
            N += 2

            # Failsafe
            if N > 200:
                raise ValueError("HalfBand.__call__: Too high order of the filter.")
