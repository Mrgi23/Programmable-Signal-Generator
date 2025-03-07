#include <cmath>
#include <stdexcept>
#include "utils.h"
#include "complexMixer.h"

using namespace std;

ComplexMixer::ComplexMixer(uint nIter, double fres) : nIter(nIter + 1), fres(fres) {
    // Compute the rotation factors 1 + j * 2 ** -i.
    factors = vector<complex<double>>(this->nIter, {0.0, 0.0});
    for (int i = 1; i < this->nIter; i++) { factors[i] = {1.0, pow(2, 1 - i)}; }

    // Add initial shift pi/2.
    factors[0] = {0.0, 1.0};
}

vector<double> ComplexMixer::operator()(
    double fshift,
    double fs,
    const vector<double>& I,
    const vector<double>& Q
) {
    if (fs <= 0) { throw invalid_argument("ComplexMixer.operator(): Sampling frequency must be positive."); }

    // Calculate the word length and the maximum phase increment.
    uint L = ceil(log2(fs / fres));
    uint Wmax = pow(2, L);

    // Compute phase increment and accumulated phases for the desired shift frequency.
    uint W = static_cast<uint>(fshift * Wmax / fs) % Wmax;
    vector<uint> Z = NCO(W, Wmax, I.size());

    // Rotate the vector.
    vector<double> Iout = CORDIC(Wmax, Z, I, Q);
    return Iout;
}

vector<double> ComplexMixer::CORDIC(
    uint Wmax,
    const vector<uint>& Z,
    const vector<double>& I,
    const vector<double>& Q
) {
    // Define the real output.
    vector<double> Iout(I.size(), 0.0);

    // Iterate through the complex input signal samples.
    for (uint i = 0; i < I.size(); i++) {
        // Define the current complex vector.
        complex<double> v = {I[i], Q[i]};

        // Define the current phase error.
        double z = static_cast<double>(Z[i]);

        // Perform CORDIC iterations.
        for (int k = 0; k < nIter; k++) {
            // Compute the current rotation angle, wrapped around 2*pi.
            double a = Wmax * atan2(factors[k].imag(), factors[k].real()) / (2 * M_PI);

            int rotation;
            // Rotate the vector forward.
            if (z > 0 && z < Wmax / 2) {
                v *= factors[k] / sqrt(1 + pow(2, -2 * k));
                rotation = -1;
            }
            // Rotate the vectore backward.
            else {
                v *= conj(factors[k]) / sqrt(1 + pow(2, -2 * k));
                rotation = 1;
            }

            // Update phase error.
            z = fmod((z + rotation * a + Wmax), Wmax);
        }
        // After finishing rotation by the desired angle, update the current real sample.
        Iout[i] = v.real();
    }
    return Iout;
}

vector<uint> ComplexMixer::NCO(uint W, uint Wmax, uint nPoints) {
    // Compute accumulated phases and wrap them around Wmax.
    vector<uint> Z = utils::linspace(0U, nPoints + 1, nPoints + 1, false);
    for (uint i = 0; i < nPoints + 1; i++) { Z[i] = (Z[i] * W) % Wmax; }
    return Z;
}
