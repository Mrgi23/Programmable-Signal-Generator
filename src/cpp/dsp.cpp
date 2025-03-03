#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "dsp.h"
#include "utils.h"

namespace dsp {
    std::vector<std::complex<double>> ComplexFFT::chirpZ(unsigned int N, const std::vector<std::complex<double>>& x) {
        // Compute convolution in next power of 2 elements atleast 2 * N - 1.
        unsigned int Nfft = 1u << static_cast<unsigned int>(std::ceil(std::log2(2 * N - 1)));

        // Define modified signal and convolution kernel.
        std::vector<std::complex<double>> xTilde(Nfft, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> k(Nfft, std::complex<double>(0.0, 0.0));

        // Compute current kernel sample and modified x sample.
        // xTilde[n] = x[n] * exp(-j x pi x n ^ 2 / N).
        // kernel[n] = exp(j x pi x n ^ 2 / N), kernel[-n] = kernel[n].
        std::complex<double> j(0.0, 1.0);
        for (int n = 0; n < N; n++) {
            double t = static_cast<double>(n * n) / static_cast<double>(N);
            k[n] = std::exp(j * M_PI * t);
            if (n > 0) { k[Nfft - n] = k[n]; }
            xTilde[n] = x[n] * std::exp(-1.0 * j * M_PI * t);
        }

        // Compute FFTs of both zero-padded sequences using the radix2 FFT.
        std::vector<std::complex<double>> Xtilde = radix2(xTilde);
        std::vector<std::complex<double>> K = radix2(k);

        // Perform convolution in frequency domain.
        std::vector<std::complex<double>> Y(Nfft, std::complex<double>(0.0, 0.0));
        for (unsigned int n = 0; n < Nfft; n++) { Y[n] = Xtilde[n] * K[n]; }

        // Perform inverse FFT (y = conj(fft(conj(Y)))) / Nfft.
        std::vector<std::complex<double>> Yconj(Nfft);
        for (unsigned int n = 0; n < Nfft; n++) { Yconj[n] = std::conj(Y[n]); }
        std::vector<std::complex<double>> y = radix2(Yconj);
        for (unsigned int n = 0; n < Nfft; n++) { y[n] = std::conj(y[n]) / static_cast<double>(Nfft); }

        // Compute chirpZ FFT of the input signal.
        // X[n] = exp(-j x pi x n ^ 2 / N) xTilde[n] * k[n].
        std::vector<std::complex<double>> X(N, std::complex<double>(0.0, 0.0));
        for (unsigned int n = 0; n < N; n++) {
            double t = static_cast<double>(n * n) / static_cast<double>(N);
            X[n] = y[n] * std::exp(-j * M_PI * t);
        }
        return X;
    }

    std::vector<std::complex<double>> ComplexFFT::dft(const std::vector<std::complex<double>>& x) {
        // Store DFT output.
        int N = x.size();
        std::vector<std::complex<double>> X(N, std::complex<double>(0.0, 0.0));

        // Perform DFT.
        // X[k] = x[n] * exp(-2 * j * pi * k * n / N).
        std::complex<double> j(0.0, 1.0);
        for (int k = 0; k < N; k++) {
            for (int n = 0; n < N; n++) {
                // Calculate twiddle factor.
                double t = static_cast<double>(k * n) / static_cast<double>(N);
                std::complex<double> W = std::exp(-2.0 * j * M_PI * t);

                // Calculate current DFT sample.
                X[k] += x[n] * W;
            }
        }
        return X;
    }

    std::vector<std::complex<double>> ComplexFFT::radix2(const std::vector<std::complex<double>>& x) {
        // Retrurn original array if it is of size 1.
        int N = x.size();
        if (N == 1) { return x; }

        // Extract even and odd elements.
        std::vector<std::complex<double>> Xeven(N / 2, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> Xodd(N / 2, std::complex<double>(0.0, 0.0));
        for (int n = 0; n < N / 2; n++) {
            Xeven[n] = x[2 * n];
            Xodd[n] = x[2 * n + 1];
        }

        // Perform radix2 recursively.
        Xeven = radix2(Xeven);
        Xodd = radix2(Xodd);

        // Store radix2 output.
        std::vector<std::complex<double>> X(N, std::complex<double>(0.0, 0.0));

        // Combine even and odd elements using radix2 algorithm.
        // X[k] = Xeven[x] + Xodd[k] * exp(-2 * j * M_PI * k / N).
        // X[k + N / 2] = Xeven[x] - Xodd[k] * exp(-2 * j * M_PI * k / N).
        std::complex<double> j(0.0, 1.0);
        for (int k = 0; k < N / 2; k++) {
            // Calculate twiddle factor.
            double t = static_cast<double>(k) / static_cast<double>(N);
            std::complex<double> W = std::exp(-2.0 * j * M_PI * t);

            // Calculate current radix2 sample.
            X[k] = Xeven[k] + W * Xodd[k];
            X[k + N / 2] = Xeven[k] - W * Xodd[k];
        }
        return X;
    }

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

    std::vector<std::complex<double>> lfilter(
        const std::vector<double>& b,
        const std::vector<std::complex<double>>& x
    ) {
        // Calculate size of the filter and of the signal.
        unsigned int nB = b.size();
        unsigned int nSignal = x.size();

        // Initialize empty filtered signal.
        std::vector<std::complex<double>> y(nSignal, std::complex<double>(0.0, 0.0));

        // Filter the input signal with the filter coefficients.
        for (unsigned int n = 0; n < nSignal; n++) {
            for (unsigned int k = 0; k < nB && k <= n; k++) { y[n] += b[k] * x[n - k]; }
        }
        return y;
    }

    std::vector<double> remez(
        unsigned int numtaps,
        const std::vector<double>& bands,
        const std::vector<double>& desired,
        const std::vector<double>& weights,
        double fs,
        unsigned int maxIter,
        unsigned int gridDensity
    ) {
        if (numtaps % 2 != 0) { throw std::invalid_argument("remez: Even number of taps required."); }

        if (bands.size() % 2 != 0) { throw std::invalid_argument("remez: Bands vector must have an even number of elements."); }
        for (unsigned int i = 0; i < bands.size() - 1; i++) {
            if (bands[i] >= bands[i + 1]) { throw std::invalid_argument("remez: Band edges must be strictly increasing."); }
        }

        unsigned int numBands = bands.size() / 2;
        if (desired.size() != numBands) { throw std::invalid_argument("remez: Desired vector must have length equal to half the number of band edges."); }

        std::vector<double> weight;
        if (weights.empty()) {
            weight = std::vector<double>(numBands, 1.0);
        }
        else {
            if (weights.size() != numBands) { throw std::invalid_argument("remez: Weight vector must have length equal to half the number of band edges."); }
            weight = weights;
        }

        if (fs <= 0.0) { throw std::invalid_argument("remez: Sampling frequency must be positive."); }

        // Normalize the band edges to [0,1] (where 1 corresponds to Nyquist = fs / 2).
        std::vector<double> normBands;
        for (double b : bands) {
            double normB = b / (fs / 2.0);
            if (normB < 0 || normB > 1) { throw std::invalid_argument("remez: Band edges must lie between 0 and fs/2."); }
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

        // For a Type I filter, set M = (numtaps - 1) / 2.
        unsigned int M = numtaps / 2;

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
