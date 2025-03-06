#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include "dac.h"
#include "dsp.h"
#include "utils.h"

using namespace std;

class TestDAC : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
        DAC converter;
};

TEST_F(TestDAC, validOutput) {
    double f = 50.0;
    double fs = 1000.0;
    vector<double> digital = utils::linspace(0.0, 1.0, static_cast<uint>(fs), 0);
    for (uint i = 0; i < digital.size(); i++) { digital[i] = cos(2 * M_PI * f * digital[i]); }
    uint nNyquist = 4;
    double Fpass = 0.4;
    double errordB = 0.025;

    vector<double> analog = converter(digital, "NRZ", nNyquist, Fpass, errordB);
    ASSERT_EQ(analog.size(), nNyquist * digital.size()) << "Invalid size of the analog signal.";
    vector<complex<double>> fftAnalog = fft(analog);
    double max = 0.0;
    for (uint i = 0; i < fftAnalog.size() / 2; i++) {
        max = abs(fftAnalog[i]) > max ? abs(fftAnalog[i]) : max;
        if (i == static_cast<uint>(fs-f) || i == static_cast<uint>(fs+f) || i == static_cast<uint>(2*fs-f)) {
            ASSERT_GT(abs(fftAnalog[i]), 0.01 * abs(fftAnalog[static_cast<uint>(f)])) << "Invalid analog signal value.";
        }
        else if (i != static_cast<uint>(f)) {
            ASSERT_LT(abs(fftAnalog[i]), 0.006 * abs(fftAnalog[static_cast<uint>(f)])) << "Invalid analog signal value.";
        }
    }
    ASSERT_EQ(abs(fftAnalog[static_cast<uint>(f)]), max) << "Invalid analog signal value.";

    analog = converter(digital, "RF", nNyquist, Fpass, errordB);
    ASSERT_EQ(analog.size(), nNyquist * digital.size()) << "Invalid size of the analog signal.";
    fftAnalog = fft(analog);
    max = 0.0;
    for (uint i = 0; i < fftAnalog.size() / 2; i++) {
        max = abs(fftAnalog[i]) > max ? abs(fftAnalog[i]) : max;
        if (i == static_cast<uint>(f) || i == static_cast<uint>(fs+f) || i == static_cast<uint>(2*fs-f)) {
            ASSERT_GT(abs(fftAnalog[i]), 0.001 * abs(fftAnalog[static_cast<uint>(fs-f)])) << "Invalid analog signal value.";
        }
        else if (i != static_cast<uint>(fs-f)) {
            ASSERT_NEAR(abs(fftAnalog[i]), 0, 1e-5) << "Invalid analog signal value.";
        }
    }
    ASSERT_EQ(abs(fftAnalog[static_cast<uint>(fs-f)]), max) << "Invalid analog signal value.";
}

TEST_F(TestDAC, invalidInput) {
    string mode = "RFF"; // Invalid reconstruction mode.
    EXPECT_THROW(converter({1}, mode), invalid_argument);

    mode = "RF";
    uint nNyquist = 5; // Indalid number of Nyquist zones for the RF mode.
    EXPECT_THROW(converter({1}, mode, nNyquist), invalid_argument);
}
