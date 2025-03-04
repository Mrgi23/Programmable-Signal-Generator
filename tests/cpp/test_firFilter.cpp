#include <gtest/gtest.h>
#include "firFilter.h"

class HalfBandTest : public ::testing::Test {
    protected:
        HalfBand halfband;
};

TEST_F(HalfBandTest, validOutput) {
    std::vector<double> coeffs = halfband(60.0, 0.2);
    ASSERT_GT(coeffs.size(), 0) << "Filter coiefficients should not be empty array.";
    ASSERT_NEAR(coeffs[coeffs.size()/2], 0.5, 1e-5) << "Center coefficient should be close to 0.5.";
    for (uint i = 0; i < coeffs.size(); i++) {
        if (i % 2 && i != coeffs.size() / 2) {
            ASSERT_EQ(coeffs[i], 0.0) << "Every odd coefficient should be 0, except center coefficient.";
        }
    }
}

TEST_F(HalfBandTest, orderIncrease) {
    std::vector<double> coeffsLower = halfband(60.0, 0.1);
    std::vector<double> coeffsHigher = halfband(60.0, 0.2);
    ASSERT_LT(coeffsLower.size(), coeffsHigher.size()) << "Higher bandpass equals higher filter order.";

    coeffsLower = halfband(60.0, 0.2);
    coeffsHigher = halfband(80.0, 0.2);
    ASSERT_LT(coeffsLower.size(), coeffsHigher.size()) << "Higher attenuation equals higher filter order.";
}

TEST_F(HalfBandTest, tooHighOrder) { EXPECT_THROW(halfband(100.0, 0.23), std::length_error); }

TEST_F(HalfBandTest, invalidInput) {
    double Fpass = -0.1; // Passband must lie between 0.0 and 0.25.
    EXPECT_THROW(halfband(60.0, Fpass), std::invalid_argument);

    Fpass = 0.3; // Passband must lie between 0.0 and 0.25.
    EXPECT_THROW(halfband(60.0, Fpass), std::invalid_argument);
}
