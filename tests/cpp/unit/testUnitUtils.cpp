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
