#ifndef COMPLEX_MIXER_H
#define COMPLEX_MIXER_H

#ifdef __cplusplus

#include <vector>
#include <complex>

class ComplexMixer {
    private:
        int nIter;
        double fres;
        std::vector<std::complex<double>> factors;

        std::vector<double> CORDIC(
            uint Wmax,
            const std::vector<uint>& Z,
            const std::vector<double>& I,
            const std::vector<double>& Q
        );
        std::vector<uint> NCO(uint W, uint Wmax, uint nPoints);
    public:
        ComplexMixer(int nIter = 13, double fres = 1.0);
        ~ComplexMixer() {}

        std::vector<double> operator()(
            double fshift,
            double fs,
            const std::vector<double>& I,
            const std::vector<double>& Q
        );
};

#endif

#endif
