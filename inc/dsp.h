#ifndef DSP_H
#define DSP_H

#ifdef __cplusplus

#include <complex>
#include <vector>

namespace dsp {
    std::vector<std::complex<double>> freqz(
        std::vector<double>& w,
        const std::vector<double>& b,
        unsigned int worN = 1024,
        double fs = 1
    );
    std::vector<double> remez(
        unsigned int numtaps,
        const std::vector<double>& bands,
        const std::vector<double>& desired,
        const std::vector<double>& weight,
        double fs = 1.0,
        unsigned int maxIter = 25,
        unsigned int gridDensity = 16
    );
}

#endif

#endif
