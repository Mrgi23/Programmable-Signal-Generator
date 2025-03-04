#include <cmath>
#include <complex>
#include <stdexcept>
#include "dsp.h"
#include "firFilter.h"
#include "utils.h"
#include <iostream>

std::vector<double> InverseSinc::operator()(double Fpass, double errordB, uint nSpec) {
    if (Fpass < 0.0 || Fpass > 0.5) { throw std::invalid_argument("InverseSinc.operator(): Passband must lie between 0.0 and 0.5."); }

    // Define spectrum of frequencies.
    std::vector<double> f = utils::linspace(0.0, Fpass, nSpec + 1);

    // Compute target frequencies & target frequency response (inverse sinc) and reduce singularity.
    std::vector<double> ftarget(2 * (nSpec + 1), 0.0);
    std::vector<double> htarget(2 * (nSpec + 1), 0.0);
    for (uint i = 0; i < nSpec + 1; i++) {
        // Target frequencies.
        ftarget[2*i] = 2 * f[i];
        ftarget[2*i+1] = 2 * (f[i] + 1e-3);

        // Target frequency response.
        htarget[2*i] = i > 0 ? (M_PI * f[i]) / sin(M_PI * f[i]) : 1.0;
        htarget[2*i+1] = i > 0 ? (M_PI * f[i]) / sin(M_PI * f[i]) : 1.0;
    }

    // Type I filter, odd number of taps.
    uint N = 1;

    // Iterate until filter is created.
    while (true) {
        // Design the filter.
        std::vector<double> b = dsp::firls(N, ftarget, htarget);

        // Compute the frequency response.
        std::vector<std::complex<double>> hinverse = dsp::freqz(f, b, nPoints);

        // Calculate sinc response and reduce singularity.
        std::vector<std::complex<double>> hsinc(nPoints, 0.0);
        for (uint i = 0; i < nPoints; i++) { hsinc[i] = i > 0 ? sin(M_PI * f[i]) / (M_PI * f[i]) : 1.0; }

        // Calculate error and error frequency response.
        double error = pow(10, std::abs(errordB) / 20);
        std::vector<double> herror;
        for (uint i = 0; i < nPoints; i++) {
            if (f[i] < Fpass) { herror.push_back(std::abs(hsinc[i] * hinverse[i])); }
        }

        // Check if the filter satisfy constrains.
        bool isCreated = true;
        for (uint i = 0; i < herror.size(); i++) {
            if (herror[i] <= 1 / error || herror[i] >= error) {
                isCreated = false;
                break;
            }
        }

        // Create the filter.
        if (isCreated) { return b; }

        // Increase the filter order.
        N += 2;

        // Failsafe.
        if (N > 200) { throw std::length_error("InverseSinc.operator(): Too high order of the filter."); }
    }
}

std::vector<double> HalfBand::operator()(double AdB, double Fpass) {
    if (Fpass < 0.0 || Fpass > 0.25) { throw std::invalid_argument("HalfBand.operator(): Passband must lie between 0.0 and 0.25."); }

    // Passband ripple.
    double deltaPass = std::pow(10.0, -std::abs(AdB) / 20.0);

    // Harris formula for initial filter order.
    uint N = static_cast<uint>(2 * std::abs(AdB) / (23 * (0.5 - 2 * Fpass)));

    // Type II filter, even number of taps.
    if (N % 2) { N += 1; }

    // Allocate memory for the frequency values.
    std::vector<double> f(nPoints, 0.0);

    // Iterate until filter is created.
    while (true) {
        // Design the filter.
        std::vector<double> b;
        try { b = dsp::remez(N, {0.0, 2 * Fpass}, {1.0}); }
        catch(const std::exception) { b = std::vector<double>(N, 0.0); }

        // Calculate the frequency response.
        std::vector<std::complex<double>> h = dsp::freqz(f, b, nPoints);

        // Calculate error and error frequency response.
        double error = 2 * deltaPass;
        std::vector<double> herror;
        for (uint i = 0; i < nPoints; i++) {
            if (f[i] < Fpass) { herror.push_back(std::abs(std::abs(h[i]) - 1.0)); }
        }

        // Check if the filter satisfy constrains.
        bool isCreated = true;
        for (uint i = 0; i < herror.size(); i++) {
            if (herror[i] > error) {
                isCreated = false;
                break;
            }
        }

        // Create the halfband filter.
        if (isCreated) {
            std::vector<double> coeffs(2 * N - 1, 0);
            for (uint i = 0; i < N; i++) { coeffs[2*i] = 0.5 * b[i]; }
            coeffs[N-1] = 0.5;
            return coeffs;
        }

        // Increase the filter order.
        N += 2;

        // Failsafe.
        if (N > 200) { throw std::length_error("HalfBand.operator(): Too high order of the filter."); }
    }
}
