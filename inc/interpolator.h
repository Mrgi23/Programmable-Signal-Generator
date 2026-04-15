#pragma once

#include "types.h"
#include "firFilter.h"
#include <complex>

class Interpolator {
    public:
        Interpolator(uint nSteps = 4, uint nPoints = 8192);
        ~Interpolator();

        uint getNSteps(void);

        std::vector<std::complex<double>> operator()(
            double AdB,
            double fmax,
            double fs,
            const std::vector<std::complex<double>>& input
        );
    private:
        std::vector<std::complex<double>> filter(
            const std::vector<double>& b,
            const std::vector<std::complex<double>>& input
        );
        std::vector<std::complex<double>> upsample(
            const std::vector<std::complex<double>>& input
        );

        uint nSteps;
        HalfBand * halfband;
};
