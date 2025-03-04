#include <gtest/gtest.h>
#include "utils.h"

TEST(UtilsTest, linspaceValidOutput) {
    uint num = 0;
    std::vector<double> vec = utils::linspace(2.0, 3.0, num);
    ASSERT_EQ(vec.size(), num) << "Invalid linear space size.";

    double start = 2.0;
    num = 1;
    vec = utils::linspace(start, 3.0, num);
    ASSERT_EQ(vec.size(), num) << "Invalid linear space size.";
    ASSERT_EQ(vec[0], start) << "Only element should be start.";

    double end = 3.0;
    num = 5;
    int endpoint = 1;
    std::vector<double> vecExpected = {2.0, 2.25, 2.5, 2.75, 3.0};

    vec = utils::linspace(start, end, num, endpoint);
    ASSERT_EQ(vec.size(), num) << "Invalid linear space size.";
    for (uint i = 0; i < num; i++) { ASSERT_EQ(vec[i], vecExpected[i]) << "Invalid slinear space value."; }

    endpoint = 0;
    vecExpected = {2.0, 2.2, 2.4, 2.6, 2.8};

    vec = utils::linspace(start, end, num, endpoint);
    ASSERT_EQ(vec.size(), num) << "Invalid linear space size.";
    for (uint i = 0; i < num; i++) { ASSERT_EQ(vec[i], vecExpected[i]) << "Invalid slinear space value."; }
}

TEST(UtilsTest, solveValidOutput) {
    std::vector<std::vector<double>> A = {{2, 1}, {5, 3}};
    std::vector<double> b = {3, 7};
    std::vector<double> xExpected = {2.0, -1.0};

    std::vector<double> x = utils::solve(A, b);
    ASSERT_EQ(x.size(), xExpected.size()) << "Invalid size of the resulting array.";
    for (uint i = 0; i < x.size(); i++) { ASSERT_NEAR(x[i], xExpected[i], 1e-5)  << "Invalid result value."; }
}

TEST(UtilsTest, solveSingularMatrix) {
    std::vector<std::vector<double>> A = {{2, 1}, {4, 2}}; // Matrix is singular or nearly singular.
    EXPECT_THROW(utils::solve(A, {7, 8}), std::runtime_error);
}

TEST(UtilsTest, solveInvalidInput) {
    std::vector<double> b = {7, 8, 5}; // Vector must have the same number of rows, as the matrix.
    EXPECT_THROW(utils::solve({{2, 1}, {1, 3}}, b), std::invalid_argument);

    std::vector<std::vector<double>> A = {{2, 1, 4}, {1, 3, 4}}; // Matrix must be square.
    EXPECT_THROW(utils::solve(A, {7, 8}), std::invalid_argument);
}
