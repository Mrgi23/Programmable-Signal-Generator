#include <algorithm>
#include <cmath>
#include <complex>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "utils.h"
#include "dsp.h"
#include "dac.h"

using namespace std;

// Define the test object.
class TestDAC : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
        DAC dac;
};

TEST_F(TestDAC, testDAC) {
    // Define the input.
    double f = 50.0;
    double fs = 1000.0;
    vector<double> digital = utils::linspace(0.0, 1.0, static_cast<uint>(fs), 0);
    for (uint i = 0; i < digital.size(); i++) { digital[i] = cos(2 * M_PI * f * digital[i]); }
    uint nNyquist = 4;
    double Fpass = 0.4;
    double errordB = 0.025;

    // Compute the result.
    vector<double> analog = dac(digital, "NRZ", nNyquist, Fpass, errordB);
    vector<complex<double>> fftAnalog = fft.fft(analog);

    // Test the result.
    vector<uint> spec = {static_cast<uint>(f)};
    for (uint i = 1; i < nNyquist / 2 + 1; i++) {
        spec.push_back(static_cast<uint>(i * fs - f));
        spec.push_back(static_cast<uint>(i * fs + f));
    }
    spec.pop_back();
    double maxValue = 0.0;
    for (uint i = 0; i < fftAnalog.size() / 2; i++) { maxValue = maxValue > abs(fftAnalog[i]) ? maxValue : abs(fftAnalog[i]); }

    ASSERT_EQ(analog.size(), nNyquist * digital.size()) << "Invalid size of the analog signal.";
    ASSERT_EQ(abs(fftAnalog[static_cast<uint>(f)]), maxValue) << "Original spectral component must be at f.";
    for (uint i = 0; i < fftAnalog.size() / 2; i++) {
        if (find(spec.begin(), spec.end(), i) != spec.end()) {
            ASSERT_GT(abs(fftAnalog[i]), 0.01 * maxValue) << "Spectral replicas must be above 1% of maximum value.";
        }
        else { ASSERT_LT(abs(fftAnalog[i]), 0.006 * maxValue) << "Noise must be below 0.6% of maximum value."; }
    }

    // Compute the result.
    analog = dac(digital, "RF", nNyquist, Fpass, errordB);
    fftAnalog = fft.fft(analog);

    // Test the result.
    maxValue = 0.0;
    for (uint i = 0; i < fftAnalog.size() / 2; i++) { maxValue = maxValue > abs(fftAnalog[i]) ? maxValue : abs(fftAnalog[i]); }

    ASSERT_EQ(analog.size(), nNyquist * digital.size()) << "Invalid size of the analog signal.";
    ASSERT_EQ(abs(fftAnalog[static_cast<uint>(fs - f)]), maxValue) << "Original spectral component must be at fs - f.";
    for (uint i = 0; i < fftAnalog.size() / 2; i++) {
        if (find(spec.begin(), spec.end(), i) != spec.end()) {
            ASSERT_GT(abs(fftAnalog[i]), 0.001 * maxValue) << "Spectral replicas must be above 0.1% of maximum value.";
        }
        else { ASSERT_NEAR(abs(fftAnalog[i]), 0, 1e-5) << "Noise must be close to 0."; }
    }
}
