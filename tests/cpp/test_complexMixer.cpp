#include <cmath>
#include <gtest/gtest.h>
#include "dsp.h"
#include "complexMixer.h"
#include "utils.h"

class TestComplexMixer : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
        ComplexMixer complexMixer;
};

TEST_F(TestComplexMixer, validOutput) {
    double fshift = 250.0;
    double fs = 1000.0;
    double f = 100.0;
    std::vector<double> I = utils::linspace(0.0, 1.0, static_cast<uint>(fs), 0);
    for (uint i = 0; i < I.size(); i++) { I[i] = cos(2 * M_PI * f * I[i]); }
    std::vector<double> Q(static_cast<uint>(fs), 0);

    std::vector<double> Iout = complexMixer(fshift, fs, I, Q);
    std::vector<std::complex<double>> fftIout = fft(Iout);
    ASSERT_EQ(fftIout.size(), I.size()) << "Invalid size of the shifted signal.";
    for (uint i = 0; i < fftIout.size() / 2; i++) {
        double sampleExpected = 0.0;
        if (i == static_cast<int>(fshift-f) || i == static_cast<int>(fshift+f)) { sampleExpected = fs / 4; }
        ASSERT_NEAR(abs(fftIout[i]), sampleExpected, 1e-5) << "Invalid shifted signal value.";
    }
}

TEST_F(TestComplexMixer, invalidInput) {
    double fs = -1000.0; // Sampling frequency must be positive.
    EXPECT_THROW(complexMixer(100.0, fs, {1}, {1}), std::invalid_argument);
}
