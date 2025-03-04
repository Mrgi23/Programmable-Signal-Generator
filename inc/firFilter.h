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

class HalfBand : public FIR {
    public:
        HalfBand(uint nPoints = 8192) : FIR(nPoints) {}
        ~HalfBand() override {}

        std::vector<double> operator()(double AdB, double Fpass);
};

#endif

#endif
