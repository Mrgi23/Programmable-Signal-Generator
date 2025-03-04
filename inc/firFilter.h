#ifndef FIR_FILTER_H
#define FIR_FILTER_H

#ifdef __cplusplus

#include <vector>

class FIR {
    protected:
        int nPoints;
    public:
        FIR(uint nPoints = 8192) : nPoints(nPoints) {}
        virtual ~FIR() {}
};

class InverseSinc : public FIR {
    public:
        InverseSinc(uint nPoints = 8192) : FIR(nPoints) {}
        ~InverseSinc() override {}

        std::vector<double> operator()(double Fpass, double errordB, uint nSpec = 16);
};

class HalfBand : public FIR {
    public:
        HalfBand(uint nPoints = 8192) : FIR(nPoints) {}
        ~HalfBand() override {}

        std::vector<double> operator()(double AdB, double Fpass);
};

#endif

#endif
