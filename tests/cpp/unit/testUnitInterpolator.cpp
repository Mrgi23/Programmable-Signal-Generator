#include <algorithm>
#include <cmath>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "dsp.h"

#define private public
#include "firFilter.h"
#include "interpolator.h"
#undef private

using namespace std;
using ::testing::_;
using ::testing::Return;

// Define the mock object.
class HalfBandMock : public HalfBand {
    public:
        HalfBandMock(int nPoints = 8192) : HalfBand(nPoints) {}
        MOCK_METHOD(vector<double>, Call, (double AdB, double Fpass), ());
        vector<double> operator()(double AdB, double Fpass) override { return Call(AdB, Fpass); }
};

// Define the test object.
class TestInterpolator : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
        Interpolator interpolator;
};

TEST_F(TestInterpolator, validOutput) {
    // Define the output.
    double fmax = 400.0;
    double fs = 1000.0;
    vector<complex<double>> input(static_cast<int>(fs), complex<double>(0.0, 0.0));
    for (uint n = 0; n < static_cast<int>(fs); n++) {
        input[n] = complex<double>(sin(2 * M_PI * fmax * static_cast<double>(n) / fs));
    }

    // Mock the HalfBand filter.
    vector<double> halfbandResult = {1.0};
    delete interpolator.halfband;
    interpolator.halfband = new HalfBandMock();
    HalfBandMock * halfbandMock = dynamic_cast<HalfBandMock*>(interpolator.halfband);
    EXPECT_CALL(*halfbandMock, Call(_, _)).Times(4).WillRepeatedly(Return(halfbandResult));

    // Compute the result.
    vector<complex<double>> output = interpolator(120.0, fmax, fs, input);
    output = fft.fft(output);
    uint N = static_cast<uint>(pow(2U, interpolator.getN()) * fs);

    // Test the result.
    vector<uint> spec;
    for (uint i  = 0; i < N / 2; i++) {
        spec.push_back(static_cast<uint>(i * fs + fmax));
        spec.push_back(static_cast<uint>((i + 1) * fs - fmax));
    }

    ASSERT_EQ(output.size(), N) << "Invalid size of the interpolated signal.";
    for (uint i = 0; i < N / 2; i++) {
        if (find(spec.begin(), spec.end(), i) != spec.end()) {
            ASSERT_NEAR(abs(output[i]), fs / 2, 1e-9) << "Spectral components must be at the multiples of fmax and fs - fmax.";
        }
        else { ASSERT_NEAR(abs(output[i]), 0, 1e-9) << "Noise must be close to 0."; }
    }
}

TEST_F(TestInterpolator, invalidInput) {
    double fs = -1000.0; // Sampling frequency must be positive.
    EXPECT_THROW(interpolator(60.0, 200.0, fs, {1}), invalid_argument);
}
