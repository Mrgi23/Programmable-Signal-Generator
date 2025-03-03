#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#ifdef __cplusplus

#include <complex>
#include <vector>
#include "firFilter.h"

class Interpolator {
    private:
        HalfBand halfband;

        std::vector<std::complex<double>> filter(
            std::vector<double> b,
            std::vector<std::complex<double>> input
        );
        std::vector<std::complex<double>> upsample(
            unsigned int n,
            std::vector<std::complex<double>> input
        );
    public:
        Interpolator() {}
        ~Interpolator() {}

        std::vector<std::complex<double>> operator()(
            double AdB,
            double fmax,
            double fs,
            std::vector<std::complex<double>> input
        );
};

#endif

#endif
