#include "dsp.h"
#include "interpolator.h"
#include <stdexcept>

using namespace std;

Interpolator::Interpolator(uint nSteps, uint nPoints) : nSteps(nSteps), halfband(new HalfBand(nPoints)) {}

Interpolator::~Interpolator()
{
    delete halfband;
}

uint Interpolator::getNSteps(void)
{
    return nSteps;
}

vector<complex<double>> Interpolator::operator()(double AdB, double fmax, double fs, const vector<complex<double>>& input)
{
    if (fs <= 0.0)
        throw invalid_argument("Interpolator.operator(): Sampling frequency must be positive.");

    vector<complex<double>> output = input;
    for (int i = 1; i < nSteps + 1; i++)
    {
        output = upsample(output);

        int factor = pow(2, i);
        double Fpass = fmax / (factor * fs);
        vector<double> b = (*halfband)(AdB, Fpass);

        output = filter(b, output);
    }
    return output; }

vector<complex<double>> Interpolator::filter(const vector<double>& b, const vector<complex<double>>& input)
{
    uint N = input.size() + b.size();
    vector<complex<double>> output(N, {0.0, 0.0});
    for (int i = 0; i < N; i++)
        output[i] = input[i % input.size()];

    output = dsp::lfilter(b, output);
    output.erase(output.begin(), output.begin() + b.size());
    return output; }

vector<complex<double>> Interpolator::upsample(const vector<complex<double>>& input)
{
    vector<complex<double>> output(2 * input.size(), {0.0, 0.0});
    for (uint i = 0; i < input.size(); i++)
        output[2 * i] = input[i];
    return output;
}
