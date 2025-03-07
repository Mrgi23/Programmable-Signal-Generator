#include <stdexcept>
#include "utils.h"
#include "dsp.h"
#include "dac.h"

using namespace std;

std::vector<double> DAC::operator()(
    const std::vector<double>& digital,
    std::string mode,
    uint nNyquist,
    double Fpass,
    double errordB
) {
    // Compute the reconstruction kernel.
    vector<double> K = kernel(mode, nNyquist);

    // Filter the input signal with the sinc compensation filter, if necessary.
    vector<double> filteredDigital;
    if (mode == "NRZ") {
        vector<double> b = inverseSinc(Fpass, errordB);
        filteredDigital = dsp::lfilter(b, digital);
    }
    else { filteredDigital = digital; }

    // Define the analog signal.
    uint N = (filteredDigital.size() - 1) * nNyquist + 1;
    vector<double> analog(N, 0.0);

    // Upsample the input signal.
    for (uint i = 0; i < N; i += nNyquist) { analog[i] = filteredDigital[i / nNyquist]; }

    // Reconstruct the digital signal as analog.
    analog = dsp::convolve(K, analog);
    return analog;
}

vector<double> DAC::kernel(std::string mode, uint nNyquist) {
    // Define kernel based on reconstruction function.
    if (mode == "NRZ") {
        // Zero-order hold.
        vector<double> K(nNyquist, 1.0);
        return K;
    }
    if (mode == "RF") {
        if (nNyquist % 2) { throw invalid_argument("DAC: Indalid number of Nyquist zones for the RF mode."); }

        // Bipolar zero-order hold.
        vector<double> K(nNyquist, 1.0);
        for (uint i = nNyquist / 2; i < nNyquist; i++) { K[i] = -1.0; }
        return K;
    }
    throw invalid_argument("DAC: Invalid reconstruction mode.");
}
