#ifndef SIGNAL_GENERATOR_H
#define SIGNAL_GENERATOR_H

#ifdef __cplusplus

#include <complex>
#include <string>
#include <vector>
#include "interpolator.h"
#include "complexMixer.h"
#include "dac.h"

class SignalGenerator {
    private:
        Interpolator * interpolator;
        ComplexMixer * complexMixer;
        DAC * dac;
    public:
        SignalGenerator(uint N = 4, uint nPoints = 8192, uint nIter = 13, double fres = 1.0);
        ~SignalGenerator();

        std::vector<double> operator()(
            const std::vector<std::complex<double>>& signal,
            double fs,
            double fmax,
            double fshift,
            std::string mode,
            double AdB = 60.0,
            uint nNyquist = 4U,
            double Fpass = 0.4,
            double errordB = 0.025
        );
};

#endif

#endif
