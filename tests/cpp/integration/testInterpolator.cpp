#include <algorithm>
#include <cmath>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "dsp.h"
#include "interpolator.h"

using namespace std;

// Define the test object.
class TestInterpolator : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
        Interpolator interpolator = Interpolator(4);
};

TEST_F(TestInterpolator, testInterpolator) {
    // Define the output.
    double fmax = 400.0;
    double fs = 1000.0;
    vector<complex<double>> input(static_cast<int>(fs), complex<double>(0.0, 0.0));
    for (uint n = 0; n < static_cast<int>(fs); n++) {
        input[n] = complex<double>(cos(2 * M_PI * fmax * static_cast<double>(n) / fs));
    }

    // Compute the result.
    vector<complex<double>> output = interpolator(120.0, fmax, fs, input);
    output = fft.fft(output);
    uint N = static_cast<uint>(pow(2U, interpolator.getNSteps()) * fs);

    // Test the result.
    ASSERT_EQ(output.size(), N) << "Invalid size of the interpolated signal.";
    for (uint i = 0; i < N / 2; i++) {
        if (i == static_cast<uint>(fmax)) {
            ASSERT_NEAR(abs(output[i]), fs / 2, 5e-3) << "Spectral component must be at fmax.";
        }
        else { ASSERT_NEAR(abs(output[i]), 0, 5e-3) << "Noise must be close to 0."; }
    }
}
