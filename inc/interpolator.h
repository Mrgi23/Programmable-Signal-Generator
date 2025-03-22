#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#ifdef __cplusplus

#include <complex>
#include <vector>
#include "firFilter.h"

class Interpolator {
    private:
        uint nSteps;
        HalfBand * halfband;

        std::vector<std::complex<double>> filter(
            const std::vector<double>& b,
            const std::vector<std::complex<double>>& input
        );
        std::vector<std::complex<double>> upsample(
            const std::vector<std::complex<double>>& input
        );
    public:
        Interpolator(uint nSteps = 4, uint nPoints = 8192) : nSteps(nSteps), halfband(new HalfBand(nPoints)) {}
        ~Interpolator() { delete halfband; }

        inline uint getNSteps(void) { return nSteps;}

        std::vector<std::complex<double>> operator()(
            double AdB,
            double fmax,
            double fs,
            const std::vector<std::complex<double>>& input
        );
};

#endif

#endif
