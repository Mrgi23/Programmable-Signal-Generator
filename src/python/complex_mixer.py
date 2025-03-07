import numpy as np

class ComplexMixer():
    def __init__(self, n_iter: int = 13, f_res: float = 1.0) -> None:
        self.n_iter = n_iter
        self.f_res = f_res

        # Compute the rotation factors 1 + j * 2 ** -i
        self.factors = 1 + 1j * 2 ** np.linspace(1, -(n_iter - 1), n_iter + 1)

        # Add initial shift pi/2.
        self.factors[0] = 1j

    def __call__(self, f_shift: float, fs: float, I_in: np.ndarray, Q_in: np.ndarray) -> np.ndarray:
        if fs <= 0:
            raise ValueError("ComplexMixer.__call__: Sampling frequency must be positive.")

        # Calculate the word length and the maximum phase increment.
        L = int(np.ceil(np.log2(fs / self.f_res)))
        W_max = 2 ** L

        # Compute phase increment and accumulated phases for the desired shift frequency.
        W = int(f_shift * W_max / fs) % W_max
        Z = self.__NCO(W, W_max, len(I_in))

        # Rotate the vector.
        I_out = self.__CORDIC(W_max, Z, I_in, Q_in)
        return I_out

    def __CORDIC(self, W_max: int, Z: np.ndarray, I_in: np.ndarray, Q_in: np.ndarray) -> np.ndarray:
        # Define the real output.
        I_out = np.zeros(len(I_in))

        # Iterate through the complex input signal samples.
        for i, (x, y, z) in enumerate(zip(I_in, Q_in, Z)):
            # Define the current complex vector.
            v = x + 1j * y

            # Perform CORDIC iterations.
            for k, factor in enumerate(self.factors):
                # Compute the current rotation angle, wrapped around 2*pi.
                a = W_max * np.angle(factor) / (2 * np.pi)

                # Rotate the vector forward.
                if z > 0 and z < W_max / 2:
                    v *= factor / np.sqrt(1 + 2 ** (-2 * k))
                    rotation = -1
                # Rotate the vector backward.
                else:
                    v *= np.conj(factor) / np.sqrt(1 + 2 ** (-2 * k))
                    rotation = 1

                # Update phase error.
                z = np.fmod(z + rotation * a + W_max, W_max)

            # After finishing rotation by the desired angle, update the current real sample.
            I_out[i] = v.real
        return I_out

    def __NCO(self, W: int, W_max: int, n_points: int) -> np.ndarray:
        # Compute accumulated phases.
        Z = W * np.arange(n_points + 1)

        # Wrap phases around 2*pi.
        Z = Z % W_max
        return Z
