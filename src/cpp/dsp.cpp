#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "dsp.h"
#include "utils.h"

namespace dsp {
    std::vector<std::complex<double>> freqz(
        std::vector<double>& w,
        const std::vector<double>& b,
        unsigned int worN,
        double fs
    ) {
        // Store frequency values.
        if (w.size() != worN) { w.resize(worN); }

        // Store frequency response.
        std::vector<std::complex<double>> h(worN);

        // Iterate through the frequency points.
        for (unsigned int i = 0; i < worN; i++) {
            double omega = M_PI * i / worN;

            // Compute the current sample of frequency value w.
            w[i] = omega * fs / (2 * M_PI);

            // Compute the current sample of frequency response H.
            for (unsigned int n = 0; n < b.size(); n++) { h[i] += b[n] * std::exp(std::complex<double>(0, -omega * n)); }
        }
        return h;
    }

    std::vector<double> remez(
        unsigned int numtaps,
        const std::vector<double>& bands,
        const std::vector<double>& desired,
        const std::vector<double>& weight,
        double fs,
        unsigned int maxIter,
        unsigned int gridDensity
    ) {
        // Ensure an even number of taps (required for a Type II design)
        if (numtaps % 2 != 0) { throw std::invalid_argument("Even number of taps required."); }
        unsigned int M = numtaps / 2;

        // Normalize the band edges to [0,1] (where 1 corresponds to Nyquist = fs / 2).
        std::vector<double> normBands;
        for (double b : bands) {
            double normB = b / (fs / 2.0);
            if (normB < 0 || normB > 1) { throw std::invalid_argument("Band edges must lie between 0 and fs/2."); }
            normBands.push_back(normB);
        }

        // Build a dense frequency grid along with desired response and weights.
        std::vector<double> grid, gridDesired, gridWeight;
        for (unsigned int i = 0; i < normBands.size() / 2; i++) {
            // Define start and stop frequency.
            double fStart = normBands[2 * i];
            double fStop = normBands[2 * i + 1];

            // Compute number of grid points in the current band.
            unsigned int nPoints = std::max(int(std::ceil(gridDensity * numtaps * (fStop - fStart))), 2);

            // Create linearly spaced values manually.
            grid = utils::linspace(fStart, fStop, nPoints);

            // Fill in the desired response and weight for this band.
            gridDesired.insert(gridDesired.end(), nPoints, desired[i]);
            gridWeight.insert(gridWeight.end(), nPoints, weight[i]);
        }

        // The number of unknown variables
        unsigned int nUnknowns = M + 1;

        // Number of extremums (at least as number of unknown libraries).
        unsigned int nExt = nUnknowns;

        // Initialize extremal frequencies by selecting nExt equally spaced points from the grid.
        std::vector<double> ext;
        {
            // Compute step in index-space.
            double step = (grid.size() - 1) / nExt;
            for (unsigned int i = 0; i < nExt; i++) {
                int idx = int(round(i * step));
                ext.push_back(grid[idx]);
            }
        }

        // Allocate vector for prototype coefficients.
        std::vector<double> bProto(M, 0.0);

        // Remez exchange iterations.
        for (unsigned int iter = 0; iter < maxIter; iter++) {
            // Initialize the interpolation matrix A and right-hand side vector b.
            std::vector<std::vector<double>> A(nExt, std::vector<double>(nUnknowns, 0.0));
            std::vector<double> b(nExt, 0.0);

            // Calculate the interpolation matrix A and right-hand side vector b.
            for (unsigned int i = 0; i < nExt; i++) {
                // Fill in cosine basis functions (scaled by 2) for each prototype coefficient.
                double f = ext[i];
                for (unsigned int k = 0; k < M; k++) { A[i][k] = 2.0 * cos(M_PI * (k + 0.5) * f); }

                // Locate the grid point closest to the current extremal frequency.
                unsigned int bestIdx = 0;
                double bestDist = fabs(grid[0] - f);
                for (unsigned int j = 1; j < grid.size(); j++) {
                    double dist = fabs(grid[j] - f);
                    if (dist < bestDist) {
                        bestDist = dist;
                        bestIdx = j;
                    }
                }
                // Set alternating error term for the Remez algorithm.
                A[i][M] = ((i % 2) == 0 ? 1.0 : -1.0) / gridWeight[bestIdx];

                // Set desired response at the extremal frequency.
                b[i] = gridDesired[bestIdx];
            }

            // Solve the linear system A * x = b.
            std::vector<double> x = utils::solve(A, b);

            // Update prototype coefficients from the solution.
            for (unsigned int i = 0; i < M; i++) { bProto[i] = x[i]; }

            // Evaluate the filter response over the entire frequency grid.
            std::vector<double> Hgrid(grid.size(), 0.0);
            for (unsigned int j = 0; j < grid.size(); j++) {
                double sum = 0.0;
                for (unsigned int k = 0; k < M; k++)
                    sum += bProto[k] * cos(M_PI * (k + 0.5) * grid[j]);
                Hgrid[j] = 2.0 * sum;
            }

            // Compute error and weighted error over the grid.
            std::vector<double> error(grid.size(), 0.0);
            std::vector<double> weightedError(grid.size(), 0.0);
            for (unsigned int j = 0; j < grid.size(); j++) {
                error[j] = gridDesired[j] - Hgrid[j];
                weightedError[j] = gridWeight[j] * error[j];
            }

            // Identify local extrema in the weighted error curve.
            std::vector<unsigned int> extIndices;
            for (unsigned int j = 1; j < grid.size() - 1; j++) {
                if ((weightedError[j] >= weightedError[j - 1] && weightedError[j] >= weightedError[j + 1]) ||
                    (weightedError[j] <= weightedError[j - 1] && weightedError[j] <= weightedError[j + 1])) {
                    extIndices.push_back(j);
                }
            }
            // Always include the endpoints.
            if (extIndices.empty() || extIndices.front() != 0) { extIndices.insert(extIndices.begin(), 0); }
            if (extIndices.back() != grid.size() - 1) { extIndices.push_back(grid.size() - 1); }

            // Select nExt grid points with the largest absolute weighted error.
            std::vector<std::pair<double, int>> errIndex;
            for (unsigned int idx : extIndices) { errIndex.push_back({fabs(weightedError[idx]), idx}); }
            sort(errIndex.begin(), errIndex.end(), [](auto& a, auto& b) { return a.first > b.first; });

            std::vector<unsigned int> newExtIndices;
            for (unsigned int i = 0; i < nExt && i < errIndex.size(); i++) { newExtIndices.push_back(errIndex[i].second); }
            sort(newExtIndices.begin(), newExtIndices.end());

            std::vector<double> newExt;
            for (int idx : newExtIndices) { newExt.push_back(grid[idx]); }

             // Check for convergence: if the extremal frequencies haven't changed significantly, exit.
            int converged = (newExt.size() == ext.size());
            for (unsigned int i = 0; i < newExt.size() && converged; i++) {
                if (fabs(newExt[i] - ext[i]) > 1e-6) { converged = false; }
            }
            ext = newExt;
            if (converged) { break; }
        }

        // Construct the full symmetric FIR filter coefficients.
        std::vector<double> bFull(numtaps, 0.0);
        for (unsigned int i = 0; i < M; i++) { bFull[i] = bProto[M - 1 - i]; }
        for (unsigned int i = 0; i < M; i++) { bFull[M + i] = bProto[i]; }
        return bFull;
    }
}
