#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "utils.h"

namespace utils {
    std::vector<double> lstsq(std::vector<std::vector<double>> A, std::vector<double> b, double lambda) {
        uint m = A.size();

        if (b.size() != m) { throw std::invalid_argument("lstsq: Vector must have the same number of rows, as the matrix."); }

        uint n = A[0].size();

        // For square A, solve linear system.
        if (m == n) { return solve(A, b); }

        // For non-square A, form the normal equations: A^T A * x = A^T b.
        std::vector<std::vector<double>> AtA(n, std::vector<double>(n, 0.0));
        for (uint i = 0; i < n; i++) {
            for (uint j = 0; j < n; j++) {
                double sum = 0.0;
                for (uint k = 0; k < m; k++) { sum += A[k][i] * A[k][j]; }
                AtA[i][j] = sum;

                // Add regularization on the diagonal to avoid singularity.
                AtA[i][j] += i == j ? lambda : 0.0;
            }
        }

        std::vector<double> Atb(n, 0.0);
        for (uint i = 0; i < n; i++) {
            double sum = 0.0;
            for (uint k = 0; k < m; k++) {
                sum += A[k][i] * b[k];
            }
            Atb[i] = sum;
        }

        // Solve the square system.
        return solve(AtA, Atb);
    }

    std::vector<double> solve(std::vector<std::vector<double>> A, std::vector<double> b) {
        uint n = A.size();

        // Check if A is square and if b has the correct size.
        if (b.size() != n) { throw std::invalid_argument("solve: Vector must have the same number of rows, as the matrix."); }
        for (const auto& row : A) {
            if (row.size() != n) { throw std::invalid_argument("solve: Matrix must be square."); }
        }

        // Forward elimination: Convert A to an upper triangular matrix.
        for (uint i = 0; i < n; i++) {
            // Partial pivoting: Find the row with the largest absolute value in the current column.
            int pivot = i;
            double maxVal = fabs(A[i][i]);
            for (uint row = i + 1; row < n; row++) {
                if (fabs(A[row][i]) > maxVal) {
                    maxVal = fabs(A[row][i]);
                    pivot = row;
                }
            }

            // If the pivot element is nearly zero, the matrix is singular.
            if (fabs(A[pivot][i]) < 1e-12) { throw std::runtime_error("solve: Matrix is singular or nearly singular."); }

            // Swap the current row with the pivot row.
            std::swap(A[i], A[pivot]);
            std::swap(b[i], b[pivot]);

            // Eliminate the entries below the pivot.
            for (uint row = i + 1; row < n; row++) {
                double factor = A[row][i] / A[i][i];
                for (uint col = i; col < n; col++)
                    A[row][col] -= factor * A[i][col];
                b[row] -= factor * b[i];
            }
        }
        // Back substitution: Solve for x in the upper triangular matrix.
        std::vector<double> x(n, 0.0);
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (uint j = i + 1; j < n; j++)
                sum += A[i][j] * x[j];
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }
}
