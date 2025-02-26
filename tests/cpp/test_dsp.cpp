#include <gtest/gtest.h>
#include "dsp.h"

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
    for (unsigned int i = 0; i < 5; i++) {
        ASSERT_NEAR(w[i], wExpected[i], 1e-5) << "Invalid frquency value.";
        ASSERT_NEAR(h[i].real(), hExpected[i].real(), 1e-5) << "Invalid frquency response, real part.";
        ASSERT_NEAR(h[i].imag(), hExpected[i].imag(), 1e-5) << "Invalid frquency response, imaginary part.";
    }
}

TEST(DSPTest, remezValidOutput) {
    std::vector<double> bands = {0.0, 0.4};
    std::vector<double> desired = {1.0};
    std::vector<double> weight = {1.0};
    std::vector<double> bExpected = {0.10753661, -0.18306401,  0.62633938,  0.62633938, -0.18306401, 0.10753661};

    std::vector<double> b = dsp::remez(6, bands, desired, weight);
    ASSERT_EQ(b.size(), bExpected.size())  << "Invalid size of the coefficients.";
    for (unsigned int i = 0; i < b.size(); i++) { ASSERT_NEAR(b[i], bExpected[i], 1e-4)  << "Invalid coefficient value."; }
}

TEST(DSPTest, remezFailedDesign) {
    std::vector<double> bands = {0.0, 0.2};
    std::vector<double> desired = {1.0};
    std::vector<double> weight = {1.0};

    EXPECT_THROW(dsp::remez(132, bands, desired, weight), std::runtime_error);
}

TEST(DSPTest, remezInvalidInput) {
    std::vector<double> bands = {0.0, 0.2};
    std::vector<double> desired = {1.0};
    std::vector<double> weight = {1.0};

    EXPECT_THROW(dsp::remez(7, bands, desired, weight), std::invalid_argument);

    bands = {0.0, 0.6};

    EXPECT_THROW(dsp::remez(6, bands, desired, weight), std::invalid_argument);
}
