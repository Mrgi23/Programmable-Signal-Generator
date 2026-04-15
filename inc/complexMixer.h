#pragma once

#include "types.h"
#include <complex>
#include <vector>

class ComplexMixer
{
    public:
        ComplexMixer(uint nIter = 13, double fres = 1.0);
        ~ComplexMixer();

        std::vector<double> operator()(double fshift, double fs, const std::vector<double>& I, const std::vector<double>& Q);
    private:
        std::vector<double> CORDIC(
            uint Wmax,
            const std::vector<double>& Z,
            const std::vector<double>& I,
            const std::vector<double>& Q
        );
        std::vector<double> NCO(double W, uint Wmax, uint nPoints);

        uint nIter;
        double fres;
        std::vector<std::complex<double>> factors;
};
