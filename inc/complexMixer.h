#pragma once

#include "types.h"
#include <complex>
#include <vector>

/**
 * @class ComplexMixer
 * @brief Performs complex frequency shifting and mixing using CORDIC and NCO algorithms.
 *
 * This class implements a configurable complex mixer used for frequency translation
 * of I/Q signals. It provides a high‑precision CORDIC rotator, a numerically controlled
 * oscillator (NCO), and a callable operator() that applies the full mixing pipeline.
 *
 * The mixer precomputes complex rotation factors based on the number of CORDIC
 * iterations and the desired frequency resolution. It supports arbitrary input
 * vector lengths and produces a real‑valued mixed output.
 */
class ComplexMixer
{
    public:
        /**
         * @brief Construct a new ComplexMixer.
         *
         * @param[in] nIter Number of CORDIC iterations (precision control).
         * @param[in] fres Frequency resolution used for NCO and factor generation.
         */
        ComplexMixer(uint nIter = 13, double fres = 1.0);

        /**
         * @brief Destructor.
         *
         * Default cleanup.
         */
        ~ComplexMixer();

        /**
         * @brief Apply complex mixing to I/Q input signals.
         *
         * Performs frequency shifting by generating a complex exponential at
         * frequency `fshift`, then mixing it with the input I/Q vectors. The
         * result is a real‑valued signal representing the frequency‑translated
         * waveform.
         *
         * @param[in] fshift Desired frequency shift.
         * @param[in] fs Sampling frequency.
         * @param[in] I In‑phase input samples.
         * @param[in] Q Quadrature input samples.
         * @return std::vector<double> Mixed output samples.
         *
         * @throws std::invalid_argument If the sampling frequency is non‑positive.
         */
        std::vector<double> operator()(double fshift, double fs, const std::vector<double>& I, const std::vector<double>& Q);

    private:
        /**
         * @brief Perform CORDIC rotation on I/Q data.
         *
         * Applies iterative CORDIC vector rotations to compute the rotated
         * signal at a given angular frequency index. Used internally by the
         * mixer to generate precise frequency shifts without trigonometric calls.
         *
         * @param[in] Wmax Maximum angular frequency index.
         * @param[in] Z Precomputed rotation factors.
         * @param[in] I In‑phase input samples.
         * @param[in] Q Quadrature input samples.
         * @return std::vector<double> Rotated output samples.
         */
        std::vector<double> CORDIC(
            uint Wmax,
            const std::vector<double>& Z,
            const std::vector<double>& I,
            const std::vector<double>& Q
        );

        /**
         * @brief Generate a numerically controlled oscillator (NCO) waveform.
         *
         * Produces a cosine‑based lookup table used for frequency mixing. The
         * table length and resolution depend on the maximum angular frequency
         * and the number of points requested.
         *
         * @param[in] W Angular frequency.
         * @param[in] Wmax Maximum angular frequency index.
         * @param[in] nPoints Number of samples to generate.
         * @return std::vector<double> NCO waveform.
         */
        std::vector<double> NCO(double W, uint Wmax, uint nPoints);

        uint nIter;
        double fres;
        std::vector<std::complex<double>> factors;
};
