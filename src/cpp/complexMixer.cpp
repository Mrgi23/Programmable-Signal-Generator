#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include "complexMixer.h"
#include "utils.h"

ComplexMixer::ComplexMixer(int nIter, double fres) : nIter(nIter + 1), fres(fres) {
    // Compute the rotation factors 1 + j * 2 ** -i.
    std::vector<double> real = utils::linspace(1.0, 1.0, this->nIter);
    std::vector<double> imag = utils::linspace(1.0, -static_cast<double>(this->nIter - 2), this->nIter);
    factors = std::vector<std::complex<double>>(this->nIter, 0.0);
    for (uint i = 1; i < this->nIter; i++) { factors[i] = {real[i], pow(2, imag[i])}; }

    // Add initial shift pi/2.
    factors[0] = std::complex<double>(0.0, 1.0);
}

std::vector<double> ComplexMixer::operator()(double fshift, double fs, std::vector<double> I, std::vector<double> Q) {
    if (fs <= 0) { throw std::invalid_argument("ComplexMixer.operator(): Sampling frequency must be positive."); }

    // Calculate the word length and the maximum phase increment.
    int L = std::ceil(std::log2(fs / fres));
    uint Wmax = pow(2, L);

    // Compute phase increment and accumulated phases for the desired shift frequency.
    uint W = static_cast<uint>(fshift * Wmax / fs) % Wmax;
    std::vector<uint> Z = NCO(W, Wmax, I.size());

    // Rotate the vector.
    std::vector<double> Iout = CORDIC(Wmax, Z, I, Q);
    return Iout;
}

std::vector<double> ComplexMixer::CORDIC(
    uint Wmax,
    std::vector<uint> Z,
    std::vector<double> I,
    std::vector<double> Q
) {
    // Store the real output.
    std::vector<double> Iout(I.size(), 0.0);

    // Iterate through the complex input signal samples.
    std::complex<double> j(0.0, 1.0);
    for (uint i = 0; i < I.size(); i++) {
        // Current complex sample.
        std::complex<double> v(I[i], Q[i]);

        // Current phase error.
        double z = static_cast<double>(Z[i]);

        // Perfomr N CORDIC iterations.
        for (int k = 0; k < nIter; k++) {
            // Compute the current rotation angle, wrapped around 2*pi.
            double a = Wmax * atan2(factors[k].imag(), factors[k].real()) / (2 * M_PI);

            int rotation;
            // Rotate the vectore forward.
            if (z > 0 && z < Wmax / 2) {
                v *= factors[k] / sqrt(1 + pow(2, -2 * k));
                rotation = -1;
            }
            // Rotate the vectore backward.
            else {
                v *= std::conj(factors[k]) / sqrt(1 + pow(2, -2 * k));;
                rotation = 1;
            }

            // Update phase error.
            z = std::fmod((z + rotation * a + Wmax), Wmax);
        }
        // After finishing rotation by the desired angle, update the current real sample.
        Iout[i] = v.real();
    }
    return Iout;
}

std::vector<uint> ComplexMixer::NCO(uint W, uint Wmax, uint nPoints) {
    // Compute accumulated phases and wrap them around Wmax.
    std::vector<uint> Z = utils::linspace(0U, (nPoints + 1), nPoints + 1, 0);
    for (uint i = 0; i < nPoints+1; i++) { Z[i] = (Z[i] * W) % Wmax; }
    return Z;
}
