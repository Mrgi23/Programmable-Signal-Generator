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
    for (unsigned int n = 0; n < static_cast<int>(fs); n++) {
        input[n] = std::complex<double>(sin(2 * M_PI * fmax * static_cast<double>(n) / fs));
    }
    unsigned int N = static_cast<unsigned int>(16 * fs);

    std::vector<std::complex<double>> output = interpolator(60.0, fmax, fs, input);
    output = fft(output, N);
    ASSERT_EQ(output.size(), N) << "Invalid size of the interpolated signal.";
    for (unsigned int i = 0; i < N / 2; i++) {
        double sampleExpected = 0.0;
        if (i == static_cast<int>(fmax)) { sampleExpected = fs / 2; }
        ASSERT_NEAR(abs(output[i]), sampleExpected, 1e-5) << "Invalid interpolated signal value.";
    }
}

TEST_F(TestInterpolator, filterValidOutput) {
    std::vector<double> b = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::complex<double>> input = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0}
    };
    std::vector<std::complex<double>> outputExpected = {
        {3.5, 0.0}, {4.5, 0.0}, {5.5, 0.0}, {6.5, 0.0}, {5.5, 0.0}, {4.5, 0.0}, {3.5, 0.0}, {2.5, 0.0}
    };

    std::vector<std::complex<double>> output = interpolator.filter(b, input);
    ASSERT_EQ(output.size(), outputExpected.size()) << "Invalid size of the filtered signal.";
    for (unsigned int i = 0; i < output.size(); i++) {
        ASSERT_NEAR(output[i].real(), outputExpected[i].real(), 1e-5) << "Invalid filtered signal value, real part.";
        ASSERT_NEAR(output[i].imag(), outputExpected[i].imag(), 1e-5) << "Invalid filtered signal value, imaginary part.";
    }
}

TEST_F(TestInterpolator, filterInvalidInput) {
    std::vector<std::complex<double>> input = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0}
    };
    EXPECT_THROW(interpolator.filter({}, input), std::invalid_argument);
}

TEST_F(TestInterpolator, upsampleValidOutput) {
    int n = 3;
    std::vector<std::complex<double>> input = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}
    };
    std::vector<std::complex<double>> outputExpected = {
        {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
        {2.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
        {3.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
        {4.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
        {5.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}
    };

    std::vector<std::complex<double>> output = interpolator.upsample(n, input);
    ASSERT_EQ(output.size(), outputExpected.size()) << "Invalid size of the upsampled signal.";
    for (int i = 0; i < output.size(); i++) {
        ASSERT_NEAR(output[i].real(), outputExpected[i].real(), 1e-5) << "Invalid upsampled signal value, real part.";
        ASSERT_NEAR(output[i].imag(), outputExpected[i].imag(), 1e-5) << "Invalid upsampled signal value, imaginary part.";
    }
}

TEST_F(TestInterpolator, upsampleInvalidInput) {
    std::vector<std::complex<double>> input = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}
    };
    EXPECT_THROW(interpolator.upsample(0, input), std::invalid_argument);
}
