#include "dsp.h"
#include <cmath>
#include <stdexcept>

using namespace std;

namespace dsp
{
    ComplexFFT::ComplexFFT() = default;

    ComplexFFT::~ComplexFFT() = default;

    vector<double> firls(
        uint numtaps,
        const vector<double>& bands,
        const vector<double>& desired,
        const vector<double>& weight,
        double fs,
        uint gridSize
    )
    {
        if (numtaps % 2 == 0)
            throw invalid_argument("firls: Odd number of taps required.");

        if (bands.size() % 2 != 0)
            throw invalid_argument("firls: Bands vector must have even length.");

        if (desired.size() != bands.size())
            throw invalid_argument("firls: Desired vector must have length equal to the number of band edges.");

        uint numBands = bands.size() / 2;
        vector<double> weights;
        if (weight.empty())
            weights = vector<double>(numBands, 1.0);
        else
        {
            if (weight.size() != numBands)
                throw invalid_argument("firls: Weight vector must have length equal to half the number of band edges.");
            weights = weight;
        }
        if (fs <= 0.0)
            throw invalid_argument("firls: Sampling frequency must be positive.");

        vector<double> normBands;
        for (double b : bands)
        {
            double normB = b / (fs / 2.0);
            if (normB < 0 || normB > 1) throw invalid_argument("firls: Band edges must lie between 0 and 1, relative to Nyquist.");
            normBands.push_back(normB);
        }

        uint M = (numtaps - 1) / 2;
        arma::vec omega = arma::linspace(0.0, M_PI, gridSize);
        arma::vec gridDesired(gridSize, arma::fill::zeros);
        arma::vec gridWeights(gridSize, arma::fill::zeros);
        for (uint i = 0; i < gridSize; i++)
        {
            double f = (omega(i) / M_PI) * (fs / 2.0);
            for (uint b = 0; b < numBands; b++)
            {
                double fstart = normBands[2 * b];
                double fstop = normBands[2 * b + 1];

                if (f >= fstart && f <= fstop)
                {
                    double d1 = desired[2 * b];
                    double d2 = desired[2 * b + 1];
                    double t = (f - fstart) / (fstop - fstart);

                    gridDesired(i) = d1 + t * (d2 - d1);
                    gridWeights(i) = weights[b];
                    break;
                }
            }
        }

        arma::mat A(gridSize, M + 1);
        for (uint j = 0; j <= M; j++)
            A.col(j) = 2 * arma::cos(omega * j);
        A.col(0).ones();

        arma::vec sqrtW = arma::sqrt(gridWeights);
        A.each_col() %= sqrtW;
        arma::vec d = gridDesired % sqrtW;
        A.diag() += 1e-12;
        arma::vec x = arma::solve(A, d);

        vector<double> b(numtaps, 0.0);
        b[M] = x(0);
        for (uint k = 1; k <= M; k++)
        {
            b[M - k] = x(k);
            b[M + k] = x(k);
        }
        return b;
    }

    vector<complex<double>> freqz(vector<double>& w, const vector<double>& b, uint worN, double fs)
    {
        if (w.size() != worN)
            w.resize(worN);

        vector<complex<double>> h(worN);
        for (uint i = 0; i < worN; i++)
        {
            double omega = M_PI * i / worN;
            w[i] = omega * fs / (2 * M_PI);
            for (uint n = 0; n < b.size(); n++)
                h[i] += b[n] * exp(complex<double>(0, -omega * n));
        }
        return h;
    }
}
