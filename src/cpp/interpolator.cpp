#include <stdexcept>
#include "dsp.h"
#include "interpolator.h"

using namespace std;

vector<complex<double>> Interpolator::operator()(
    double AdB,
    double fmax,
    double fs,
    const vector<complex<double>>& input
) {
    if (fs <= 0.0) { throw invalid_argument("Interpolator.operator(): Sampling frequency must be positive."); }

    // Initialize the output signal.
    vector<complex<double>> output = input;

    // Propagate output signal through the interpolation 4 times.
    for (int i = 0; i < 4; i++) {
        // Upsample previous output signal by factor 2.
        output = upsample(2, output);

        // Calculate FIR filter coefficients.
        int factor = pow(2, i);
        double Fpass = fmax / (factor * fs);
        vector<double> b = halfband(AdB, Fpass);

        // Filter the upsampled signal.
        output = filter(b, output);
    }
    return output;
}

vector<complex<double>> Interpolator::filter(const vector<double>& b, const vector<complex<double>>& input) {
    // Prepare signal for filtering by expanding its size for len(b) elements.
    vector<complex<double>> output(input.size() + b.size(), complex<double>(0.0, 0.0));
    for (int i = 0; i < input.size() + b.size(); i++) { output[i] = input[i % input.size()]; }

    // Filter signal.
    output = dsp::lfilter(b, output);

    // Remove additional elements of the signal.
    output.erase(output.begin(), output.begin() + b.size());
    return output;
}

vector<complex<double>> Interpolator::upsample(uint n, const vector<complex<double>>& input) {
    // Initialize empty output with the valid size.
    vector<complex<double>> output(n * input.size(), complex<double>(0.0, 0.0));

    // Upsample signal by adding n-1 zeros between every element.
    for (uint i = 0; i < input.size(); i++) { output[i * n] = input[i]; }
    return output;
}
