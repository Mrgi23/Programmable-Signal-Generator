#include <cmath>
#include <complex>
#include <stdexcept>
#include "dsp.h"
#include "firFilter.h"

std::vector<double> HalfBand::operator()(double AdB, double Fpass) {
    // Specifications are not valid.
    if (Fpass < 0.0 || Fpass > 0.25) { throw std::invalid_argument("Passband must lie between 0.0 and 0.25."); }

    // Passband ripple.
    double deltaPass = std::pow(10.0, -std::abs(AdB) / 20.0);

    // Harris formula for initial filter order.
    unsigned int N = static_cast<unsigned int>(2 * std::abs(AdB) / (23 * (0.5 - 2 * Fpass)));

    // Ensure N is even (liquidDSP requires even length for Kaiser filter)
    // Force even number of taps, since scipy.remez uses numtaps instead of order.
    if (N % 2) { N += 1; }

    // Allocate memory for the frequency values.
    std::vector<double> w(nPoints, 0.0);

    // Iterate until filter is created.
    while (1) {
        // Design the filter.
        std::vector<double> b;
        try { b = dsp::remez(N, {0.0, 2*Fpass}, {1.0}); }
        catch(const std::exception) { b = std::vector<double>(N, 0.0); }

        // Calculate the frequency response & the amplitude.
        std::vector<std::complex<double>> h = dsp::freqz(w, b, nPoints);

        // Check if the filter satisfy constrains.
        int isCreated = 1;
        for (unsigned int i = 0; i < nPoints; i++) {
            if (w[i] < Fpass && std::abs(std::abs(h[i]) - 1) > 2 * deltaPass) {
                isCreated = 0;
                break;
            }
        }

        // Create the halfband filter
        if (isCreated) {
            std::vector<double> coeffs(2*N-1, 0);
            for (unsigned int i = 0; i < N; i++) { coeffs[2*i] = 0.5 * b[i]; }
            coeffs[N-1] = 0.5;
            return coeffs;
        }

        // Increase the filter order.
        N += 2;

        // Failsafe
        if (N > 200) { throw std::length_error("Too high order of the filter."); }
    }
}
