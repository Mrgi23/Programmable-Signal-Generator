#include "signalGenerator.h"

using namespace std;

SignalGenerator::SignalGenerator(
    uint nSteps,
    uint nPoints,
    uint nIter,
    double fres
) : interpolator(new Interpolator(nSteps, nPoints)), complexMixer(new ComplexMixer(nIter, fres)), dac(new DAC(nPoints)) {}

SignalGenerator::~SignalGenerator() {
    delete interpolator;
    delete complexMixer;
    delete dac;
}

vector<double> SignalGenerator::operator()(
    const std::vector<std::complex<double>>& signal,
    double fs,
    double fmax,
    double fshift,
    std::string mode,
    double AdB,
    uint nNyquist,
    double Fpass,
    double errordB
) {
    // Interpolate the input signal.
    vector<complex<double>> interpolated = (*interpolator)(AdB, fmax, fs, signal);
    uint scale = pow(2U, interpolator->getNSteps());

    // Shift the interpolated signal.
    vector<double> I(interpolated.size(), 0.0);
    vector<double> Q(interpolated.size(), 0.0);
    for (uint i  = 0; i < interpolated.size(); i++) {
        I[i] = interpolated[i].real();
        Q[i] = interpolated[i].imag();
    }
    vector<double> shifted = (*complexMixer)(fshift, scale * fs, I, Q);

    // Convert the digital signal to analog.
    vector<double> analog = (*dac)(shifted, mode, nNyquist, Fpass, errordB);
    return analog;
}
