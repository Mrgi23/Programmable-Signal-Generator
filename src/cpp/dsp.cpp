#include <cmath>
#include <stdexcept>
#include "dsp.h"

using namespace std;

namespace dsp {
    vector<double> firls(
        uint numtaps,
        const vector<double>& bands,
        const vector<double>& desired,
        const vector<double>& weight,
        double fs,
        uint gridSize
    ) {
        if (numtaps % 2 == 0) { throw invalid_argument("firls: Odd number of taps required."); }

        if (bands.size() % 2 != 0) { throw invalid_argument("firls: Bands vector must have even length."); }

        if (desired.size() != bands.size()) { throw invalid_argument("firls: Desired vector must have length equal to the number of band edges."); }

        uint numBands = bands.size() / 2;
        vector<double> weights;
        if (weight.empty()) {
            weights = vector<double>(numBands, 1.0);
        }
        else {
            if (weight.size() != numBands) {
                throw invalid_argument("firls: Weight vector must have length equal to half the number of band edges.");
            }
            weights = weight;
        }

        if (fs <= 0.0) { throw invalid_argument("firls: Sampling frequency must be positive."); }

        // Normalize the band edges to [0, 1] (where 1 corresponds to Nyquist = fs / 2).
        vector<double> normBands;
        for (double b : bands) {
            double normB = b / (fs / 2.0);
            if (normB < 0 || normB > 1) { throw invalid_argument("firls: Band edges must lie between 0 and 1, relative to Nyquist."); }
            normBands.push_back(normB);
        }

        // For a Type I filter, set M = (numtaps - 1) / 2.
        uint M = (numtaps - 1) / 2;

        // Define a frequency grid. Sample from 0 to pi (radians).
        arma::vec omega = arma::linspace(0.0, M_PI, gridSize);

        // Compute a dense frequency grid along with desired response and weights.
        arma::vec gridDesired(gridSize, arma::fill::zeros);
        arma::vec gridWeights(gridSize, arma::fill::zeros);
        for (uint i = 0; i < gridSize; i++) {
            // For each grid point, map the radian frequency w to Hz.
            double f = (omega(i) / M_PI) * (fs / 2.0);

            // Loop over each band (each defined by two consecutive elements in bands/desired).
            for (uint b = 0; b < numBands; b++) {
                double fstart = normBands[2 * b];
                double fstop = normBands[2 * b + 1];

                // Linear interpolation between the band edge desired values.
                if (f >= fstart && f <= fstop) {
                    double d1 = desired[2 * b];
                    double d2 = desired[2 * b + 1];
                    double t = (f - fstart) / (fstop - fstart);

                    // Frequencies outside any specified band get zero desired response and zero weight.
                    gridDesired(i) = d1 + t * (d2 - d1);
                    gridWeights(i) = weights[b];
                    break;
                }
            }
        }

        // Compute the design matrix A.
        // For a symmetric FIR filter, the frequency response (ignoring linear phase delay) is:
        // F(w) ≈ a0 + 2 * sum_{k=1}^{M} a[k] cos(k w)
        // Each row i of A is: [1, 2*cos(w[i]), 2*cos(2w[i]), ..., 2*cos(M w[i])].
        arma::mat A(gridSize, M + 1);
        for (uint j = 0; j <= M; j++) {
            A.col(j) = 2 * arma::cos(omega * j);
        }
        A.col(0).ones();

        // Apply the weights: multiply each row of A and gridDesired by sqrt(weight)
        arma::vec sqrtW = arma::sqrt(gridWeights);
        A.each_col() %= sqrtW;
        arma::vec d = gridDesired % sqrtW;

        // Add the regularization parameter to avoid singular matrix.
        A.diag() += 1e-12;

        // Solve the weighted least-squares problem: A * x = d.
        arma::vec x = arma::solve(A, d);

        // Compute the full symmetric FIR filter coefficients.
        vector<double> b(numtaps, 0.0);
        b[M] = x(0);
        for (uint k = 1; k <= M; k++) {
            b[M - k] = x(k);
            b[M + k] = x(k);
        }
        return b;
    }

    vector<complex<double>> freqz(vector<double>& w, const vector<double>& b, uint worN, double fs) {
        // Define the frequency values.
        if (w.size() != worN) { w.resize(worN); }

        // Define frequency response.
        vector<complex<double>> h(worN);

        // Iterate through the frequency points.
        for (uint i = 0; i < worN; i++) {
            double omega = M_PI * i / worN;

            // Compute the current sample of frequency value w.
            w[i] = omega * fs / (2 * M_PI);

            // Compute the current sample of frequency response H.
            for (uint n = 0; n < b.size(); n++) { h[i] += b[n] * exp(complex<double>(0, -omega * n)); }
        }
        return h;
    }
}
