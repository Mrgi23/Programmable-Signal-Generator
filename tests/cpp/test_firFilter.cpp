#include <gtest/gtest.h>
#include "firFilter.h"

class FIRTest : public ::testing::Test {
    protected:
        HalfBand halfband;
        InverseSinc inverseSinc;
};

TEST_F(FIRTest, inverseSincValidOutput) {
    std::vector<double> b = inverseSinc(0.4, 0.025);
    for (uint i = 0; i < b.size() / 2; i++) { ASSERT_EQ(b[i], b[b.size()-1-i]) << "Filter must be symetric."; }
}

TEST_F(FIRTest, inverseSincOrderIncrease) {
    std::vector<double> coeffsLower = inverseSinc(0.4, 0.025);
    std::vector<double> coeffsHigher = inverseSinc(0.45, 0.025);
    ASSERT_LT(coeffsLower.size(), coeffsHigher.size()) << "Higher bandpass equals higher filter order.";

    coeffsLower = inverseSinc(0.4, 0.25);
    coeffsHigher = inverseSinc(0.4, 0.025);
    ASSERT_LT(coeffsLower.size(), coeffsHigher.size()) << "Lower error equals higher filter order.";
}

TEST_F(FIRTest, inverseSincTooHighOrder) { EXPECT_THROW(inverseSinc(0.49, 0.0025), std::length_error); }

TEST_F(FIRTest, inverseSincInvalidInput) {
    double Fpass = -0.1; // Passband must lie between 0.0 and 0.5.
    EXPECT_THROW(halfband(60.0, Fpass), std::invalid_argument);

    Fpass = 0.6; // Passband must lie between 0.0 and 0.5.
    EXPECT_THROW(halfband(60.0, Fpass), std::invalid_argument);
}

TEST_F(FIRTest, halfBandValidOutput) {
    std::vector<double> coeffs = halfband(60.0, 0.2);
    ASSERT_NEAR(coeffs[coeffs.size()/2], 0.5, 1e-5) << "Center coefficient should be close to 0.5.";
    for (uint i = 0; i < coeffs.size(); i++) {
        if (i % 2 && i != coeffs.size() / 2) {
            ASSERT_EQ(coeffs[i], 0.0) << "Every odd coefficient should be 0, except center coefficient.";
        }
        if (i < coeffs.size() / 2) {
            ASSERT_EQ(coeffs[i], coeffs[coeffs.size()-1-i]) << "Filter must be symetric.";
        }
    }
}

TEST_F(FIRTest, halfBandOrderIncrease) {
    std::vector<double> coeffsLower = halfband(60.0, 0.1);
    std::vector<double> coeffsHigher = halfband(60.0, 0.2);
    ASSERT_LT(coeffsLower.size(), coeffsHigher.size()) << "Higher bandpass equals higher filter order.";

    coeffsLower = halfband(60.0, 0.2);
    coeffsHigher = halfband(80.0, 0.2);
    ASSERT_LT(coeffsLower.size(), coeffsHigher.size()) << "Higher attenuation equals higher filter order.";
}

TEST_F(FIRTest, halfBandTooHighOrder) { EXPECT_THROW(halfband(200.0, 0.23), std::length_error); }

TEST_F(FIRTest, halfBandInvalidInput) {
    double Fpass = -0.1; // Passband must lie between 0.0 and 0.25.
    EXPECT_THROW(halfband(60.0, Fpass), std::invalid_argument);

    Fpass = 0.3; // Passband must lie between 0.0 and 0.25.
    EXPECT_THROW(halfband(60.0, Fpass), std::invalid_argument);
}
