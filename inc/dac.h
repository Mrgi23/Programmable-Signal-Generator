#ifndef DAC_H
#define DAC_H

#ifdef __cplusplus

#include <string>
#include <vector>
#include "firFilter.h"

class DAC {
    private:
        InverseSinc inverseSinc;

        std::vector<double> kernel(std::string mode, uint nNyquist);
    public:
        DAC() {}
        ~DAC() {}

        std::vector<double> operator()(
            const std::vector<double>& digital,
            std::string mode,
            uint nNyquist = 4U,
            double Fpass = 0.4,
            double errordB = 0.025
        );
};

#endif

#endif
