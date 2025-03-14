#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include "utils.h"
#include "dsp.h"
#include "complexMixer.h"

using namespace std;

// Define the test object.
class TestComplexMixer : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
        ComplexMixer complexMixer;
};

TEST_F(TestComplexMixer, validOutput) {
    // Define the input.
    double fshift = 150.0;
    double fs = 1000.0;
    double f = 300.0;
    vector<double> I = utils::linspace(0.0, 1.0, static_cast<uint>(fs), false);
    for (uint i = 0; i < I.size(); i++) { I[i] = cos(2 * M_PI * f * I[i]); }
    vector<double> Q(static_cast<uint>(fs), 0);

    // Compute the result.
    vector<double> Iout = complexMixer(fshift, fs, I, Q);
    vector<complex<double>> fftIout = fft.fft(Iout);

    // Test the result.
    vector<uint> spec = {
        static_cast<uint>(f - fshift),
        static_cast<uint>(f + fshift)
    };

    ASSERT_EQ(fftIout.size(), I.size()) << "Invalid size of the shifted signal.";
    for (uint i = 0; i < fftIout.size() / 2; i++) {
        if (find(spec.begin(), spec.end(), i) != spec.end()) {
            ASSERT_NEAR(abs(fftIout[i]), fs / 4, 5e-2) << "Spectral components must be at f - fshift and f + fshift.";
        }
        else { ASSERT_NEAR(abs(fftIout[i]), 0, 5e-2) << "Noise must be close to 0."; }
    }
}

TEST_F(TestComplexMixer, invalidInput) {
    double fs = -1000.0; // Sampling frequency must be positive.
    EXPECT_THROW(complexMixer(200.0, fs, {1}, {1}), invalid_argument);
}
