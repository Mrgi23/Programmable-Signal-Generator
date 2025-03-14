#include <gtest/gtest.h>
#include "dsp.h"

using namespace std;

// Define the test object.
class TestFFT : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
};

TEST_F(TestFFT, fftValidOutput) {
    {
        // Define the input and the expected result.
        uint N = 3;
        vector<complex<double>> x = {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}};
        vector<complex<double>> XExpected = {{6.0, 0.0}, {-1.5, 0.8660254}, {-1.5, -0.8660254}};

        // Compute the result.
        vector<complex<double>> X = fft.fft(x, N);

        // Test the result.
        ASSERT_EQ(X.size(), XExpected.size()) << "Invalid size of the FFT signal.";
        for (uint i = 0; i < X.size(); i++) {
            ASSERT_NEAR(X[i].real(), XExpected[i].real(), 1e-8) << "Invalid FFT signal value, real part.";
            ASSERT_NEAR(X[i].imag(), XExpected[i].imag(), 1e-8) << "Invalid FFT signal value, imaginary part.";
        }
    }

    {
        // Define the input and the expected result.
        vector<complex<double>> x = {1.0, 2.0, 3.0, 4.0};
        vector<complex<double>> XExpected = {{10.0, 0.0}, {-2.0, 2.0}, {-2.0, 0.0}, {-2.0, -2.0}};

        // Compute the result.
        vector<complex<double>> X = fft.fft(x);

        // Test the result.
        ASSERT_EQ(X.size(), XExpected.size()) << "Invalid size of the FFT signal.";
        for (uint i = 0; i < X.size(); i++) {
            ASSERT_EQ(X[i].real(), XExpected[i].real()) << "Invalid FFT signal value, real part.";
            ASSERT_EQ(X[i].imag(), XExpected[i].imag()) << "Invalid FFT signal value, imaginary part.";
        }
    }

}

TEST_F(TestFFT, fftshiftValidOutput) {
    {
        // Define the input and the expected result.
        vector<complex<double>> x = {{0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}};
        vector<complex<double>> shiftedExpected = {{3.0, 0.0}, {4.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}};

        // Compute the result.
        vector<complex<double>> shifted = fft.fftshift(x);

        // Test the result.
        ASSERT_EQ(shifted.size(), shiftedExpected.size()) << "Invalid size of the FFT shifted signal.";
        for (uint i = 0; i < shifted.size(); i++) {
            ASSERT_EQ(shifted[i].real(), shiftedExpected[i].real()) << "Invalid FFT shifted signal value, real part.";
            ASSERT_EQ(shifted[i].imag(), shiftedExpected[i].imag()) << "Invalid FFT shifted signal value, imaginary part.";
        }
    }

    {
        // Define the input and the expected result.
        vector<uint> x = {0U, 1U, 2U, 3U, 4U, 5U};
        vector<uint> shiftedExpected = {3U, 4U, 5U, 0U, 1U, 2U};

        // Compute the result.
        vector<uint> shifted = fft.fftshift(x);

        // Test the result.
        ASSERT_EQ(shifted.size(), shiftedExpected.size()) << "Invalid size of the FFT shifted signal.";
        for (uint i = 0; i < shifted.size(); i++) {
            ASSERT_EQ(shifted[i], shiftedExpected[i]) << "Invalid FFT shifted signal value.";
        }
    }
}

TEST(TestDSP, convolveValidOutput) {
    // Define the input and the expected result.
    vector<double> h = {1.0, 1.0, 0.0, 0.0};
    vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    vector<double> yExpected = {1.0, 3.0, 5.0, 7.0, 9.0, 5.0, 0.0, 0.0};

    // Compute the result.
    vector<double> y = dsp::convolve(h, x);

    // Test the result.
    ASSERT_EQ(y.size(), yExpected.size())  << "Invalid size of the convoluted signal.";
    for (uint i = 0; i < y.size(); i++) { ASSERT_EQ(y[i], yExpected[i])  << "Invalid convoluted signal value."; }
}

TEST(TestDSP, firlsValidOutput) {
    // Define the input and the expected result.
    int numtaps = 5;
    vector<double> bands = {0.0, 0.2, 0.4, 0.5};
    vector<double> desired = {1.0, 0.75, 0.0, 0.0};
    vector<double> bExpected = {0.219047675, 0.077019991, 0.39487779, 0.077019991, 0.219047675};

    // Compute the result.
    vector<double> b = dsp::firls(numtaps, bands, desired);

    // Test the result.
    ASSERT_EQ(b.size(), bExpected.size())  << "Invalid size of the coefficients.";
    for (uint i = 0; i < b.size(); i++) { ASSERT_NEAR(b[i], bExpected[i], 1e-4)  << "Invalid coefficient value."; }
}

TEST(TestDSP, firlsInvalidInput) {
    uint numtaps = 6; // Odd number of taps required.
    EXPECT_THROW(dsp::firls(numtaps, {1}, {1}), invalid_argument);

    vector<double> bands = {0.0, 0.25, 0.5}; // Bands vector must have even length.
    EXPECT_THROW(dsp::firls(5, bands, {1}), invalid_argument);

    vector<double> desired = {1.0}; // Desired vector must have length equal to the number of band edges.
    EXPECT_THROW(dsp::firls(5, {0.0, 0.2, 0.3, 0.5}, desired), invalid_argument);

    vector<double> weights = {1.0}; // Weight vector must have length equal to half the number of band edges.
    EXPECT_THROW(dsp::firls(5, {0.0, 0.2, 0.3, 0.5}, {1.0, 1.0, 0.0, 0.0}, weights), invalid_argument);

    double fs = -1000.0; // Sampling frequency must be positive.
    EXPECT_THROW(dsp::firls(5, {0.0, 0.2, 0.3, 0.5}, {1.0, 1.0, 0.0, 0.0}, {1.0, 0.0}, fs), invalid_argument);

    bands = {0.0, 0.2, 0.5, 1.1}; // Band edges must lie between 0 and 1, relative to Nyquist.
    EXPECT_THROW(dsp::firls(5, bands, {1.0, 0.0}, {1.0, 0.0}), invalid_argument);
}

TEST(TestDSP, freqzValidOutput) {
    // Define the input and the expected result.
    vector<double> w(0, 0.0);
    vector<double> b = {0.1, 0.2, 0.3};
    vector<double> wExpected = {0.0, 0.1, 0.2, 0.3, 0.4};
    vector<complex<double>> hExpected = {
        {0.6, 0}, {0.3545085, -0.40287401}, {-0.0809017, -0.36654688}, {-0.2045085, -0.01387573}, {0.0309017, 0.1677599}
    };

    // Compute the result.
    vector<complex<double>> h = dsp::freqz(w, b, 5);

    // Test the result.
    ASSERT_EQ(w.size(), wExpected.size()) << "Invalid size of the frquencies.";
    ASSERT_EQ(h.size(), hExpected.size()) << "Invalid size of the frquency response.";
    for (uint i = 0; i < 5; i++) {
        ASSERT_NEAR(w[i], wExpected[i], 1e-5) << "Invalid frquency value.";
        ASSERT_NEAR(h[i].real(), hExpected[i].real(), 1e-8) << "Invalid frquency response, real part.";
        ASSERT_NEAR(h[i].imag(), hExpected[i].imag(), 1e-8) << "Invalid frquency response, imaginary part.";
    }
}

TEST(TestDSP, lfilterValidOutput) {
    // Define the input and the expected result.
    vector<double> b = {0.25, 0.25, 0.25, 0.25};
    vector<complex<double>> x = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0}
    };
    vector<complex<double>> yExpected = {
        {0.25, 0.0}, {0.75, 0.0}, {1.5, 0.0}, {2.5, 0.0}, {3.5, 0.0}, {4.5, 0.0}, {5.5, 0.0}, {6.5, 0.0}
    };

    // Compute the result.
    vector<complex<double>> y = dsp::lfilter(b, x);

    // Test the result.
    ASSERT_EQ(y.size(), yExpected.size())  << "Invalid size of the filtered signal.";
    for (uint i = 0; i < y.size(); i++) {
        ASSERT_EQ(y[i].real(), yExpected[i].real()) << "Invalid filtered signal value, real part.";
        ASSERT_EQ(y[i].imag(), yExpected[i].imag()) << "Invalid filtered signal value, imaginary part.";
    }
}
