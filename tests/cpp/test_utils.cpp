#include <gtest/gtest.h>
#include "utils.h"

using namespace std;

TEST(TestUtils, linspaceValidOutput) {
    // Define the input.
    uint num = 0;

    // Compute the result.
    vector<double> vec = utils::linspace(2.0, 3.0, num);

    // Test the result.
    ASSERT_EQ(vec.size(), num) << "Invalid linear space size.";

    // Define the input.
    double start = 2.0;
    num = 1;

    // Compute the result.
    vec = utils::linspace(start, 3.0, num);

    // Test the result.
    ASSERT_EQ(vec.size(), num) << "Invalid size of the linear space.";
    ASSERT_EQ(vec[0], start) << "Only element should be start.";

    // Define the input and the expected result.
    double end = 3.0;
    num = 5;
    bool endpoint = true;
    vector<double> vecExpected = {2.0, 2.25, 2.5, 2.75, 3.0};

    // Compute the result.
    vec = utils::linspace(start, end, num, endpoint);

    // Test the result.
    ASSERT_EQ(vec.size(), num) << "Invalid size of the linear space.";
    for (uint i = 0; i < num; i++) { ASSERT_EQ(vec[i], vecExpected[i]) << "Invalid linear space value."; }

    // Define the input and the expected result.
    endpoint = false;
    vecExpected = {2.0, 2.2, 2.4, 2.6, 2.8};

    // Compute the result.
    vec = utils::linspace(start, end, num, endpoint);

    // Test the result.
    ASSERT_EQ(vec.size(), num) << "Invalid size of the linear space.";
    for (uint i = 0; i < num; i++) { ASSERT_EQ(vec[i], vecExpected[i]) << "Invalid linear space value."; }
}

TEST(TestUtils, lstsqValidOutput) {
    // Define the input and the expected result.
    vector<vector<double>> A = {{2, 1, 4}, {1, 3, 4}};
    vector<double> b = {7, 8};
    vector<double> xExpected = {0.46666667, 0.73333333, 1.33333333};

    // Compute the result.
    vector<double> x = utils::lstsq(A, b);

    // Test the result.
    ASSERT_EQ(x.size(), xExpected.size()) << "Invalid size of the lstsq array.";
    for (uint i = 0; i < x.size(); i++) { ASSERT_NEAR(x[i], xExpected[i], 1e-5)  << "Invalid lstsq value."; }
}

TEST(TestUtils, lstsqInvalidInput) {
    vector<double> b = {7, 8, 5}; // Vector must have the same number of rows, as the matrix.
    EXPECT_THROW(utils::lstsq({{2, 1}, {1, 3}}, b), invalid_argument);
}

TEST(TestUtils, solveValidOutput) {
    // Define the input and the expected result.
    vector<vector<double>> A = {{2, 1}, {5, 3}};
    vector<double> b = {3, 7};
    vector<double> xExpected = {2.0, -1.0};

    // Compute the result.
    vector<double> x = utils::solve(A, b);

    // Test the result.
    ASSERT_EQ(x.size(), xExpected.size()) << "Invalid size of the solve array.";
    for (uint i = 0; i < x.size(); i++) { ASSERT_NEAR(x[i], xExpected[i], 1e-5)  << "Invalid solve value."; }
}

TEST(TestUtils, solveSingularMatrix) {
    vector<vector<double>> A = {{2, 1}, {4, 2}}; // Matrix is singular or nearly singular.
    EXPECT_THROW(utils::solve(A, {7, 8}), runtime_error);
}

TEST(TestUtils, solveInvalidInput) {
    vector<double> b = {7, 8, 5}; // Vector must have the same number of rows, as the matrix.
    EXPECT_THROW(utils::solve({{2, 1}, {1, 3}}, b), invalid_argument);

    vector<vector<double>> A = {{2, 1, 4}, {1, 3, 4}}; // Matrix must be square.
    EXPECT_THROW(utils::solve(A, {7, 8}), invalid_argument);
}
