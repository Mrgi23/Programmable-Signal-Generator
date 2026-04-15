#include "utils.h"
#include "complexMixer.h"
#include <cmath>
#include <stdexcept>

using namespace std;

ComplexMixer::ComplexMixer(uint nIter, double fres) : nIter(nIter + 1), fres(fres)
{
    factors = vector<complex<double>>(this->nIter, {0.0, 0.0});
    for (int i = 1; i < this->nIter; i++) { factors[i] = {1.0, pow(2, 1 - i)}; }
    factors[0] = {0.0, 1.0};
}

ComplexMixer::~ComplexMixer() = default;

vector<double> ComplexMixer::operator()(double fshift, double fs, const vector<double>& I, const vector<double>& Q)
{
    if (fs <= 0) { throw invalid_argument("ComplexMixer.operator(): Sampling frequency must be positive."); }

    uint L = ceil(log2(fs / fres));
    uint Wmax = pow(2, L);

    double W = fmod(fshift * Wmax / fs, Wmax);
    vector<double> Z = NCO(W, Wmax, I.size());

    vector<double> Iout = CORDIC(Wmax, Z, I, Q);
    return Iout;
}

vector<double> ComplexMixer::CORDIC(uint Wmax, const vector<double>& Z, const vector<double>& I, const vector<double>& Q)
{
    vector<double> Iout(I.size(), 0.0);
    for (uint i = 0; i < I.size(); i++)
    {
        complex<double> v = {I[i], Q[i]};
        double z = Z[i];

        for (int k = 0; k < nIter; k++)
        {
            double a = Wmax * atan2(factors[k].imag(), factors[k].real()) / (2 * M_PI);

            int rotation;
            if (z > 0 && z < Wmax / 2)
            {
                v *= factors[k] / sqrt(1 + pow(2, -2 * k));
                rotation = -1;
            }
            else
            {
                v *= conj(factors[k]) / sqrt(1 + pow(2, -2 * k));
                rotation = 1;
            }

            z = fmod((z + rotation * a + Wmax), Wmax);
        }
        Iout[i] = v.real();
    }

    return Iout;
}

vector<double> ComplexMixer::NCO(double W, uint Wmax, uint nPoints)
{
    vector<double> Z = utils::linspace(0.0, static_cast<double>(nPoints + 1), nPoints + 1, false);
    for (uint i = 0; i < nPoints + 1; i++) { Z[i] = fmod(Z[i] * W, Wmax); }
    return Z;
}
