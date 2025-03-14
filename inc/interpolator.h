#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#ifdef __cplusplus

#include <complex>
#include <vector>
#include "firFilter.h"

class Interpolator {
    private:
        uint N;
        HalfBand * halfband;

        std::vector<std::complex<double>> filter(
            const std::vector<double>& b,
            const std::vector<std::complex<double>>& input
        );
        std::vector<std::complex<double>> upsample(
            uint n,
            const std::vector<std::complex<double>>& input
        );
    public:
        Interpolator(uint N = 4, uint nPoints = 8192) : N(N), halfband(new HalfBand(nPoints)) {}
        ~Interpolator() { delete halfband; }

        inline uint getN(void) { return N;}

        std::vector<std::complex<double>> operator()(
            double AdB,
            double fmax,
            double fs,
            const std::vector<std::complex<double>>& input
        );
};

#endif

#endif
