#pragma once

#include "types.h"
#include <vector>

class FIR
{
    public:
        FIR(uint nPoints = 8192);
        virtual ~FIR();
    protected:
        int nPoints;
};

class InverseSinc : public FIR
{
    public:
        InverseSinc(uint nPoints = 8192);
        ~InverseSinc() override;

        virtual std::vector<double> operator()(double Fpass, double errordB, uint nSpec = 16);
};

class HalfBand : public FIR
{
    public:
        HalfBand(uint nPoints = 8192);
        ~HalfBand() override;

        virtual std::vector<double> operator()(double AdB, double Fpass);
};
