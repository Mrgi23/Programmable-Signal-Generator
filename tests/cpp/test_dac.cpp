#include <algorithm>
#include <cmath>
#include <complex>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "utils.h"
#include "dsp.h"

#define private public
#include "firFilter.h"
#include "dac.h"
#undef private

using namespace std;
using ::testing::_;
using ::testing::Return;

class InverseSincMock : public InverseSinc {
    public:
        InverseSincMock(int nPoints = 8192) : InverseSinc(nPoints) {}
        MOCK_METHOD(vector<double>, Call, (double AdB, double Fpass, uint nSpec), ());
        vector<double> operator()(double AdB, double Fpass, uint nSpec = 16) override { return Call(AdB, Fpass, nSpec); }
};

class TestDAC : public ::testing::Test {
    protected:
        dsp::ComplexFFT fft;
        DAC dac;
};

TEST_F(TestDAC, validOutput) {
    // Define the input.
    double f = 50.0;
    double fs = 1000.0;
    vector<double> digital = utils::linspace(0.0, 1.0, static_cast<uint>(fs), 0);
    for (uint i = 0; i < digital.size(); i++) { digital[i] = cos(2 * M_PI * f * digital[i]); }
    uint nNyquist = 4;
    double Fpass = 0.4;
    double errordB = 0.025;

    // Mock the InverseSinc filter.
    vector<double> inverseSincResult = {1.0};
    delete dac.inverseSinc;
    dac.inverseSinc = new InverseSincMock();
    InverseSincMock * inverseSincMock = dynamic_cast<InverseSincMock*>(dac.inverseSinc);
    EXPECT_CALL(*inverseSincMock, Call(_, _, _)).WillOnce(Return(inverseSincResult));

    // Compute the result.
    vector<double> analog = dac(digital, "NRZ", nNyquist, Fpass, errordB);
    vector<complex<double>> fftAnalog = fft(analog);

    // Test the result.
    vector<uint> spec = {static_cast<uint>(f)};
    for (uint i = 1; i < 5; i++) {
        spec.push_back(static_cast<uint>(i * fs - f));
        spec.push_back(static_cast<uint>(i * fs + f));
    }
    spec.pop_back();
    ASSERT_EQ(analog.size(), nNyquist * digital.size()) << "Invalid size of the analog signal.";
    double max = 0.0;
    for (uint i = 0; i < fftAnalog.size() / 2; i++) {
        max = abs(fftAnalog[i]) > max ? abs(fftAnalog[i]) : max;
        if (find(spec.begin(), spec.end(), i) != spec.end()) {
            ASSERT_GT(abs(fftAnalog[i]), 0.01 * abs(fftAnalog[static_cast<uint>(f)])) << "Spectral replicas must be above 1% of maximum value.";
        }
        else {
            ASSERT_LT(abs(fftAnalog[i]), 0.006 * abs(fftAnalog[static_cast<uint>(f)])) << "Noise must be below 0.6% of maximum value.";
        }
    }
    ASSERT_EQ(abs(fftAnalog[static_cast<uint>(f)]), max) << "Maximum value must be at the frqeuncy f.";

    // Compute the result.
    analog = dac(digital, "RF", nNyquist, Fpass, errordB);
    fftAnalog = fft(analog);

    // Test the result.
    ASSERT_EQ(analog.size(), nNyquist * digital.size()) << "Invalid size of the analog signal.";
    max = 0.0;
    for (uint i = 0; i < fftAnalog.size(); i++) {
        max = abs(fftAnalog[i]) > max ? abs(fftAnalog[i]) : max;
        if (find(spec.begin(), spec.end(), i) != spec.end()) {
            ASSERT_GT(abs(fftAnalog[i]), 0.01 * abs(fftAnalog[static_cast<uint>(f)])) << "Spectral replicas must be above 1% of maximum value.";
        }
        else { ASSERT_NEAR(abs(fftAnalog[i]), 0, 1e-5) << "Noise must be close to 0."; }
    }
    ASSERT_EQ(abs(fftAnalog[static_cast<uint>(fs-f)]), max) << "Maximum value must be at the frqeuncy fs-f.";
}

TEST_F(TestDAC, invalidInput) {
    string mode = "RFF"; // Invalid reconstruction mode.
    EXPECT_THROW(dac({1}, mode), invalid_argument);

    mode = "RF";
    uint nNyquist = 5; // Indalid number of Nyquist zones for the RF mode.
    EXPECT_THROW(dac({1}, mode, nNyquist), invalid_argument);
}
