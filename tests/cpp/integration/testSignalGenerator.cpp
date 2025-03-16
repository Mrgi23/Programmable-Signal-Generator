#include <cmath>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "utils.h"
#include "dsp.h"

#define private public
#include "signalGenerator.h"
#undef private

using namespace std;

// Define the test object.
class TestSignalGenerator : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
        SignalGenerator signalGenerator;
};

TEST_F(TestSignalGenerator, testSignalGenerator) {
    // Define the input.
    uint M = pow(2U, signalGenerator.interpolator->getNSteps());
    double fs = 60.0;
    double fmin = 12.0;
    double fmax = 24.0;
    double fshiftNRZ = 3 * fs * M / 8;
    double fshiftRF = 6 * fs * M / 8;
    uint nNyquist = 4;

    vector<complex<double>> signal;
    utils::readFile("../../data/testSignal.txt", signal);
    uint N = signal.size();

    // Compute the result.
    vector<double> analog = signalGenerator(signal, fs, fmax, fshiftNRZ, "NRZ", 60.0, nNyquist);
    vector<complex<double>> fftAnalog = fft.fft(analog);

    // Test the result.
    vector<uint> specMax;
    vector<uint> spec;
    double f0 = fshiftNRZ + fmin;
    double f1 = fshiftNRZ + fmax;
    double maxValue = 0.0;
    for (uint i = 0; i < fftAnalog.size() / 2; i++) {
        maxValue = maxValue > abs(fftAnalog[i]) ? maxValue : abs(fftAnalog[i]);
        if (i >= static_cast<uint>(f0 / fs * N) && i < static_cast<uint>(f1 / fs * N)) { specMax.push_back(i); }
        else if (
            (i > static_cast<uint>((M * fs - f1) / fs * N) && i <= static_cast<uint>((M * fs - f0) / fs * N)) ||
            (i >= static_cast<uint>((M * fs + f0) / fs * N) && i < static_cast<uint>((M * fs + f1) / fs * N)) ||
            (i > static_cast<uint>((2 * M * fs - f1) / fs * N) && i <= static_cast<uint>((2 * M * fs - f0) / fs * N))
        ) { spec.push_back(i); }
    }

    ASSERT_EQ(analog.size(), nNyquist * M * N) << "Invalid size of the analog signal.";
    for (uint i = 0; i < fftAnalog.size() / 2; i++) {
        if (find(specMax.begin(), specMax.end(), i) != specMax.end()) {
            ASSERT_NEAR(abs(fftAnalog[i]), maxValue, 5e-2) << "Original spectral component must be between f_shift + f_min and f_shift + f_max.";
        }
        else if (find(spec.begin(), spec.end(), i) != spec.end()) {
            ASSERT_GT(abs(fftAnalog[i]), 0.01 * maxValue) << "Spectral replicas must be above 1% of maximum value.";
        }
        else { ASSERT_LT(abs(fftAnalog[i]), 0.006 * maxValue) << "Noise must be below 0.6% of maximum value."; }
    }

    // Compute the result.
    analog = signalGenerator(signal, fs, fmax, fshiftRF, "RF", 60.0, nNyquist);
    fftAnalog = fft.fft(analog);

    // Test the result.
    specMax.clear();
    spec.clear();
    f0 = fshiftRF + fmin;
    f1 = fshiftRF + fmax;
    maxValue = 0.0;
    for (uint i = 0; i < fftAnalog.size() / 2; i++) {
        maxValue = maxValue > abs(fftAnalog[i]) ? maxValue : abs(fftAnalog[i]);
        if (i >= static_cast<uint>(f0 / fs * N) && i < static_cast<uint>(f1 / fs * N)) { specMax.push_back(i); }
        else if (
            (i > static_cast<uint>((M * fs - f1) / fs * N) && i <= static_cast<uint>((M * fs - f0) / fs * N)) ||
            (i >= static_cast<uint>((M * fs + f0) / fs * N) && i < static_cast<uint>((M * fs + f1) / fs * N)) ||
            (i > static_cast<uint>((2 * M * fs - f1) / fs * N) && i <= static_cast<uint>((2 * M * fs - f0) / fs * N))
        ) { spec.push_back(i); }
    }

    ASSERT_EQ(analog.size(), nNyquist * M * N) << "Invalid size of the analog signal.";
    for (uint i = 0; i < fftAnalog.size() / 2; i++) {
        if (find(specMax.begin(), specMax.end(), i) != specMax.end()) {
            ASSERT_NEAR(abs(fftAnalog[i]), maxValue, 5e-2) << "Original spectral component must be between f_shift + f_min and f_shift + f_max.";
        }
        else if (find(spec.begin(), spec.end(), i) != spec.end()) {
            ASSERT_GT(abs(fftAnalog[i]), 0.001 * maxValue) << "Spectral replicas must be above 0.1% of maximum value.";
        }
        else { ASSERT_NEAR(abs(fftAnalog[i]), 0, 1e-3) << "Noise must be close to 0."; }

    }
}
