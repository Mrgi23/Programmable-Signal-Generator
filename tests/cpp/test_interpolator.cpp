#include <cmath>
#include <gtest/gtest.h>
#include "dsp.h"
#include "interpolator.h"

class TestInterpolator : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
        Interpolator interpolator;
};

TEST_F(TestInterpolator, validOutput) {
    double fmax = 200.0;
    double fs = 1000.0;
    std::vector<std::complex<double>> input(static_cast<int>(fs), std::complex<double>(0.0, 0.0));
    for (uint n = 0; n < static_cast<int>(fs); n++) {
        input[n] = std::complex<double>(sin(2 * M_PI * fmax * static_cast<double>(n) / fs));
    }
    uint N = static_cast<uint>(16 * fs);

    std::vector<std::complex<double>> output = interpolator(60.0, fmax, fs, input);
    output = fft(output, N);
    ASSERT_EQ(output.size(), N) << "Invalid size of the interpolated signal.";
    for (uint i = 0; i < N / 2; i++) {
        double sampleExpected = 0.0;
        if (i == static_cast<int>(fmax)) { sampleExpected = fs / 2; }
        ASSERT_NEAR(abs(output[i]), sampleExpected, 1e-5) << "Invalid interpolated signal value.";
    }
}

TEST_F(TestInterpolator, invalidInput) {
    double fs = -1000.0; // Sampling frequency must be positive.
    EXPECT_THROW(interpolator(60.0, 200, fs, {1}), std::invalid_argument);
}