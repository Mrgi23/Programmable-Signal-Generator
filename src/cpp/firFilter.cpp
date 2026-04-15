#include "utils.h"
#include "dsp.h"
#include "firFilter.h"
#include <cmath>
#include <complex>
#include <liquid/liquid.h>
#include <stdexcept>

using namespace std;

FIR::FIR(uint nPoints) : nPoints(nPoints) {}

FIR::~FIR() = default;

InverseSinc::InverseSinc(uint nPoints) : FIR(nPoints) {}

InverseSinc::~InverseSinc() = default;

vector<double> InverseSinc::operator()(double Fpass, double errordB, uint nSpec)
{
    if (Fpass <= 0.0 || Fpass >= 0.5)
        throw invalid_argument("InverseSinc.operator(): Passband must lie between 0.0 and 0.5.");

    vector<double> f = utils::linspace(0.0, Fpass, nSpec + 1);

    vector<double> fTarget(2 * (nSpec + 1), 0.0);
    vector<double> hTarget(2 * (nSpec + 1), 0.0);
    for (uint i = 0; i < nSpec + 1; i++)
    {
        fTarget[2*i] = 2 * f[i];
        fTarget[2*i+1] = 2 * (f[i] + 1e-3);

        hTarget[2*i] = i > 0 ? (M_PI * f[i]) / sin(M_PI * f[i]) : 1.0;
        hTarget[2*i+1] = i > 0 ? (M_PI * f[i]) / sin(M_PI * f[i]) : 1.0;
    }

    uint N = 1;
    while (true)
    {
        vector<double> b = dsp::firls(N, fTarget, hTarget);
        vector<complex<double>> hinverse = dsp::freqz(f, b, nPoints);
        double error = pow(10, abs(errordB) / 20);

        bool isCreated = true;
        for (uint i = 0; i < nPoints; i++)
            if (f[i] < Fpass)
            {
                double hsinc = i > 0 ? sin(M_PI * f[i]) / (M_PI * f[i]) : 1.0;
                if(abs(hsinc * hinverse[i]) < 1 / error || abs(hsinc * hinverse[i]) > error)
                {
                    isCreated = false;
                    break;
                }
            }
        if (isCreated)
            return b;

        N += 2;
        if (N > 200)
            throw length_error("InverseSinc.operator(): Too high order of the filter.");
    }
}

HalfBand::HalfBand(uint nPoints) : FIR(nPoints) {}

HalfBand::~HalfBand() = default;

vector<double> HalfBand::operator()(double AdB, double Fpass)
{
    if (Fpass <= 0.0 || Fpass >= 0.25)
        throw invalid_argument("HalfBand.operator(): Passband must lie between 0.0 and 0.25.");

    double deltaPass = pow(10.0, -abs(AdB) / 20.0);
    uint N = static_cast<uint>(abs(AdB) / (46 * (0.5 - 2 * Fpass)));
    if (N % 2)
        N += 1;

    vector<double> f(nPoints, 0.0);
    while (true)
    {
        vector<float> bands = {0.0f, static_cast<float>(2 * Fpass)};
        vector<float> des = {1.0f, 0.0f};
        vector<float> weights = {1.0};
        vector<liquid_firdespm_wtype> wtype = {LIQUID_FIRDESPM_FLATWEIGHT};
        vector<float> c(N, 0.0f);

        int firOK = firdespm_run(
            N, 1, bands.data(), des.data(), weights.data(), wtype.data(), LIQUID_FIRDESPM_BANDPASS, c.data()
        );

        if (firOK == LIQUID_OK && !isnan(c[1]))
        {
            vector<double> b(c.begin(), c.end());
            vector<complex<double>> h = dsp::freqz(f, b, nPoints);
            double error = 2 * deltaPass;

            bool isCreated = true;
            for (uint i = 0; i < nPoints; i++)
                if (f[i] < 2 * Fpass && abs(abs(h[i]) - 1.0) > error)
                {
                    isCreated = false;
                    break;
                }

            // Compute the halfband filter.
            if (isCreated)
            {
                vector<double> coeffs(2 * N - 1, 0.0);
                for (uint i = 0; i < N; i++) { coeffs[2 * i] = 0.5 * b[i]; }
                coeffs[N - 1] = 0.5;
                return coeffs;
            }
        }

        N += 2;
        if (N > 200)
            throw length_error("HalfBand.operator(): Too high order of the filter.");
    }
}
