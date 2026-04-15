#include "utils.h"
#include "dsp.h"
#include "dac.h"
#include <stdexcept>

using namespace std;

DAC::DAC(uint nPoints) : inverseSinc(new InverseSinc(nPoints)) {}

DAC::~DAC() { delete inverseSinc; }

vector<double> DAC::operator()(const vector<double>& digital, string mode, uint nNyquist, double Fpass, double errordB)
{
    vector<double> K = kernel(mode, nNyquist);

    vector<double> filteredDigital;
    if (mode == "NRZ")
    {
        vector<double> b = (*inverseSinc)(Fpass, errordB);
        filteredDigital = dsp::lfilter(b, digital);
    }
    else
        filteredDigital = digital;

    uint N = (filteredDigital.size() - 1) * nNyquist + 1;
    vector<double> analog(N, 0.0);

    for (uint i = 0; i < N; i += nNyquist)
        analog[i] = filteredDigital[i / nNyquist];

    analog = dsp::convolve(K, analog);
    return analog;
}

vector<double> DAC::kernel(string mode, uint nNyquist)
{
    if (mode == "NRZ")
    {
        vector<double> K(nNyquist, 1.0);
        return K;
    }

    if (mode == "RF")
    {
        if (nNyquist % 2)
            throw invalid_argument("DAC.kernel: Invalid number of Nyquist zones for the RF mode.");

        vector<double> K(nNyquist, 1.0);
        for (uint i = nNyquist / 2; i < nNyquist; i++)
            K[i] = -1.0;
        return K;
    }

    throw invalid_argument("DAC.kernel: Invalid reconstruction mode.");
}
