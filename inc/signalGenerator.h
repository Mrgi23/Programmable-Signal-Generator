#pragma once

#include "types.h"
#include "interpolator.h"
#include "complexMixer.h"
#include "dac.h"
#include <string>

class SignalGenerator {
    public:
        SignalGenerator(uint nSteps = 4, uint nPoints = 8192, uint nIter = 13, double fres = 1.0);
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
     private:
        Interpolator * interpolator;
        ComplexMixer * complexMixer;
        DAC * dac;
};
