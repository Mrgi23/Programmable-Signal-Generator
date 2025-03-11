#include <gtest/gtest.h>
#include "firFilter.h"

using namespace std;

class TestFIR : public ::testing::Test {
    protected:
        HalfBand halfband;
        InverseSinc inverseSinc;
};

TEST_F(TestFIR, inverseSincValidOutput) {
    // Compute the result.
    vector<double> b = inverseSinc(0.4, 0.025);

    // Test the result.
    for (uint i = 0; i < b.size() / 2; i++) { ASSERT_EQ(b[i], b[b.size()-1-i]) << "Filter must be symetric."; }
}

TEST_F(TestFIR, inverseSincOrderIncrease) {
    // Compute the result.
    vector<double> bLower = inverseSinc(0.4, 0.025);
    vector<double> bHigher = inverseSinc(0.45, 0.025);

    // Test the result.
    ASSERT_LT(bLower.size(), bHigher.size()) << "Higher bandpass equals higher filter order.";

    // Compute the result.
    bLower = inverseSinc(0.4, 0.25);
    bHigher = inverseSinc(0.4, 0.025);

    // Test the result.
    ASSERT_LT(bLower.size(), bHigher.size()) << "Lower error equals higher filter order.";
}

TEST_F(TestFIR, inverseSincTooHighOrder) {
    // Too high order of the filter.
    EXPECT_THROW(inverseSinc(0.49, 0.0025), length_error);
}

TEST_F(TestFIR, inverseSincInvalidInput) {
    double Fpass = -0.1; // Passband must lie between 0.0 and 0.5.
    EXPECT_THROW(halfband(60.0, Fpass), invalid_argument);

    Fpass = 0.6; // Passband must lie between 0.0 and 0.5.
    EXPECT_THROW(halfband(60.0, Fpass), invalid_argument);
}

TEST_F(TestFIR, halfBandValidOutput) {
    // Compute the result.
    vector<double> coeffs = halfband(120.0, 0.2);

    // Test the result.
    ASSERT_NEAR(coeffs[coeffs.size()/2], 0.5, 5e-3) << "Center coefficient should be close to 0.5.";
    for (uint i = 0; i < coeffs.size(); i++) {
        if (i % 2 && i != coeffs.size() / 2) {
            ASSERT_EQ(coeffs[i], 0.0) << "Every odd coefficient should be 0, except center coefficient.";
        }
        if (i < coeffs.size() / 2) {
            ASSERT_EQ(coeffs[i], coeffs[coeffs.size()-1-i]) << "Filter must be symetric.";
        }
    }
}

TEST_F(TestFIR, halfBandOrderIncrease) {
    // Compute the result.
    vector<double> coeffsLower = halfband(120.0, 0.1);
    vector<double> coeffsHigher = halfband(120.0, 0.2);

    // Test the result.
    ASSERT_LT(coeffsLower.size(), coeffsHigher.size()) << "Higher bandpass equals higher filter order.";

    // Compute the result.
    coeffsLower = halfband(120.0, 0.2);
    coeffsHigher = halfband(130.0, 0.2);

    // Test the result.
    ASSERT_LT(coeffsLower.size(), coeffsHigher.size()) << "Higher attenuation equals higher filter order.";
}

TEST_F(TestFIR, halfBandTooHighOrder) {
    // Too high order of the filter.
    EXPECT_THROW(halfband(270.0, 0.24), length_error);
}

TEST_F(TestFIR, halfBandInvalidInput) {
    double Fpass = -0.1; // Passband must lie between 0.0 and 0.25.
    EXPECT_THROW(halfband(120.0, Fpass), invalid_argument);

    Fpass = 0.3; // Passband must lie between 0.0 and 0.25.
    EXPECT_THROW(halfband(120.0, Fpass), invalid_argument);
}
