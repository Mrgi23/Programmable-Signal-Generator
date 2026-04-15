#pragma once

#include "types.h"
#include "firFilter.h"
#include <string>

class DAC
{
    public:
        DAC(uint nPoints = 8192);
        ~DAC();

        std::vector<double> operator()(
            const std::vector<double>& digital,
            std::string mode,
            uint nNyquist = 4U,
            double Fpass = 0.4,
            double errordB = 0.025
        );
    private:
        std::vector<double> kernel(std::string mode, uint nNyquist);

        InverseSinc * inverseSinc;
};
