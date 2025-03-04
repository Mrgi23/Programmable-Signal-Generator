#include <gtest/gtest.h>
#include "dsp.h"

class DSPTestFFT : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
};

TEST_F(DSPTestFFT, validOutput) {
    uint N = 3;
    std::vector<std::complex<double>> x = {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}};
    std::vector<std::complex<double>> XExpected = {{6.0, 0.0}, {-1.5, 0.8660254}, {-1.5, -0.8660254}};

    std::vector<std::complex<double>> X = fft(x, N);
    ASSERT_EQ(X.size(), XExpected.size()) << "Invalid size of the FFT signal.";
    for (uint i = 0; i < X.size(); i++) {
        ASSERT_NEAR(X[i].real(), XExpected[i].real(), 1e-5) << "Invalid FFT signal value, real part.";
        ASSERT_NEAR(X[i].imag(), XExpected[i].imag(), 1e-5) << "Invalid FFT signal value, imaginary part.";
    }

    XExpected = {{10.0, 0.0}, {-2.0, 2.0}, {-2.0, 0.0}, {-2.0, -2.0}};

    X = fft(x);
    ASSERT_EQ(X.size(), XExpected.size()) << "Invalid size of the FFT signal.";
    for (uint i = 0; i < X.size(); i++) {
        ASSERT_NEAR(X[i].real(), XExpected[i].real(), 1e-5) << "Invalid FFT signal value, real part.";
        ASSERT_NEAR(X[i].imag(), XExpected[i].imag(), 1e-5) << "Invalid FFT signal value, imaginary part.";
    }

    x = {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}};
    XExpected = {{15.0, 0.0}, {-2.5, 3.4409548}, {-2.5, 0.81229924}, {-2.5, -0.81229924}, {-2.5, -3.4409548}};

    X = fft(x);
    ASSERT_EQ(X.size(), XExpected.size()) << "Invalid size of the FFT signal.";
    for (uint i = 0; i < X.size(); i++) {
        ASSERT_NEAR(X[i].real(), XExpected[i].real(), 1e-5) << "Invalid FFT signal value, real part.";
        ASSERT_NEAR(X[i].imag(), XExpected[i].imag(), 1e-5) << "Invalid FFT signal value, imaginary part.";
    }
}

TEST(DSPTest, freqzValidOutput) {
    std::vector<double> w(0, 0.0);
    std::vector<double> b = {0.1, 0.2, 0.3};
    std::vector<double> wExpected = {0.0, 0.1, 0.2, 0.3, 0.4};
    std::vector<std::complex<double>> hExpected = {
        {0.6, 0}, {0.3545085, -0.40287401}, {-0.0809017, -0.36654688}, {-0.2045085, -0.01387573}, {0.0309017, 0.1677599}
    };

    std::vector<std::complex<double>> h = dsp::freqz(w, b, 5);
    ASSERT_EQ(w.size(), wExpected.size()) << "Invalid size of the frquencies.";
    ASSERT_EQ(h.size(), hExpected.size()) << "Invalid size of the frquency response.";
    for (uint i = 0; i < 5; i++) {
        ASSERT_NEAR(w[i], wExpected[i], 1e-5) << "Invalid frquency value.";
        ASSERT_NEAR(h[i].real(), hExpected[i].real(), 1e-5) << "Invalid frquency response, real part.";
        ASSERT_NEAR(h[i].imag(), hExpected[i].imag(), 1e-5) << "Invalid frquency response, imaginary part.";
    }
}

TEST(DSPTest, lfilterValidOutput) {
    std::vector<double> b = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::complex<double>> x = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0}
    };
    std::vector<std::complex<double>> yExpected = {
        {0.25, 0.0}, {0.75, 0.0}, {1.5, 0.0}, {2.5, 0.0}, {3.5, 0.0}, {4.5, 0.0}, {5.5, 0.0}, {6.5, 0.0}
    };

    std::vector<std::complex<double>> y = dsp::lfilter(b, x);
    ASSERT_EQ(y.size(), yExpected.size())  << "Invalid size of the filtered signal.";
    for (uint i = 0; i < y.size(); i++) {
        ASSERT_NEAR(y[i].real(), yExpected[i].real(), 1e-5) << "Invalid filtered signal value, real part.";
        ASSERT_NEAR(y[i].imag(), yExpected[i].imag(), 1e-5) << "Invalid filtered signal value, imaginary part.";
    }
}

TEST(DSPTest, remezValidOutput) {
    std::vector<double> bands = {0.0, 0.4};
    std::vector<double> desired = {1.0};
    std::vector<double> weight = {1.0};
    std::vector<double> bExpected = {0.10753661, -0.18306401,  0.62633938,  0.62633938, -0.18306401, 0.10753661};

    std::vector<double> b = dsp::remez(6, bands, desired, weight);
    ASSERT_EQ(b.size(), bExpected.size())  << "Invalid size of the coefficients.";
    for (uint i = 0; i < b.size(); i++) { ASSERT_NEAR(b[i], bExpected[i], 1e-4)  << "Invalid coefficient value."; }
}

TEST(DSPTest, remezInvalidInput) {
    uint numtaps = 7; // Even number of taps required.
    EXPECT_THROW(dsp::remez(numtaps, {1}, {1}), std::invalid_argument);

    std::vector<double> bands = {0.0, 0.25, 0.5}; // Bands vector must have an even number of elements.
    EXPECT_THROW(dsp::remez(6, bands, {1}), std::invalid_argument);

    bands = {0.0, 0.2, 0.1, 0.5}; // Band edges must be strictly increasing.
    EXPECT_THROW(dsp::remez(6, bands, {1}), std::invalid_argument);

    std::vector<double> desired = {1.0}; // Desired vector must have length equal to half the number of band edges.
    EXPECT_THROW(dsp::remez(6, {0.0, 0.2, 0.3, 0.5}, desired), std::invalid_argument);

    std::vector<double> weights = {1.0}; // Weight vector must have length equal to half the number of band edges.
    EXPECT_THROW(dsp::remez(6, {0.0, 0.2, 0.3, 0.5}, {1.0, 0.0}, weights), std::invalid_argument);

    double fs = -1000.0; // Sampling frequency must be positive.
    EXPECT_THROW(dsp::remez(6, {0.0, 0.2, 0.3, 0.5}, {1.0, 0.0}, {1.0, 0.0}, fs), std::invalid_argument);

    bands = {0.0, 0.2, 0.5, 0.6}; // Band edges must lie between 0 and fs/2.
    EXPECT_THROW(dsp::remez(6, bands, {1.0, 0.0}, {1.0, 0.0}), std::invalid_argument);
}
