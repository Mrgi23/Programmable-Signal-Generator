#include <stdexcept>
#include "dsp.h"
#include "interpolator.h"

std::vector<std::complex<double>> Interpolator::operator()(
    double AdB,
    double fmax,
    double fs,
    std::vector<std::complex<double>> input
) {
    if (fs <= 0.0) { throw std::invalid_argument("Interpolator.operator(): Sampling frequency must be positive."); }

    // Initialize the output signal.
    std::vector<std::complex<double>> output = input;

    // Propagate output signal through the interpolation 4 times.
    for (int i = 0; i < 4; i++) {
        // Upsample previous output signal by factor 2.
        output = upsample(2, output);

        // Calculate FIR filter coefficients.
        int factor = pow(2, i);
        double Fpass = fmax / (factor * fs);
        std::vector<double> b = halfband(AdB, Fpass);

        // Filter the upsampled signal.
        output = filter(b, output);
    }
    return output;
}

std::vector<std::complex<double>> Interpolator::filter(
    std::vector<double> b,
    std::vector<std::complex<double>> input
) {
    // Prepare signal for filtering by expanding its size for len(b) elements.
    std::vector<std::complex<double>> output(input.size() + b.size(), std::complex<double>(0.0, 0.0));
    for (int i = 0; i < input.size() + b.size(); i++) { output[i] = input[i % input.size()]; }

    // Filter signal.
    output = dsp::lfilter(b, output);

    // Remove additional elements of the signal.
    output.erase(output.begin(), output.begin() + b.size());
    return output;
}

std::vector<std::complex<double>> Interpolator::upsample(
    unsigned int n,
    std::vector<std::complex<double>> input
) {
    // Initialize empty output with the valid size.
    std::vector<std::complex<double>> output(n * input.size(), std::complex<double>(0.0, 0.0));

    // Upsample signal by adding n-1 zeros between every element.
    for (unsigned int i = 0; i < input.size(); i++) { output[i * n] = input[i]; }
    return output;
}
