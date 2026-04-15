#pragma once

#include "types.h"
#include "firFilter.h"
#include <complex>

/**
 * @class Interpolator
 * @brief Multistage complex-signal interpolator using halfband filters.
 *
 * This class performs multistage interpolation of complex baseband signals.
 * Each stage upsamples the signal by a factor of two and applies a halfband
 * low‑pass filter to suppress spectral images. The number of interpolation
 * stages is configurable, and the halfband filter is designed adaptively
 * based on the desired attenuation and passband edge.
 *
 * The interpolator is suitable for sample‑rate conversion, digital upconversion,
 * and multirate DSP pipelines where efficient halfband filtering is preferred.
 */
class Interpolator {
    public:
        /**
         * @brief Construct a new Interpolator.
         *
         * @param[in] nSteps Number of interpolation stages (each ×2).
         * @param[in] nPoints Number of frequency‑domain evaluation points for filter design.
         */
        Interpolator(uint nSteps = 4, uint nPoints = 8192);

        /**
         * @brief Destructor.
         *
         * Frees the internal halfband filter designer.
         */
        ~Interpolator();

        /**
         * @brief Get the number of interpolation stages.
         * @return uint Number of ×2 interpolation steps.
         */
        uint getNSteps(void);

        /**
         * @brief Perform multistage interpolation on a complex input signal.
         *
         * The signal is repeatedly upsampled by a factor of two and filtered
         * using a halfband low‑pass filter. The passband edge is adjusted at
         * each stage based on the maximum signal frequency `fmax` and the
         * sampling frequency `fs`. The halfband filter is designed using the
         * specified attenuation `AdB`.
         *
         * The sampling frequency must be strictly positive. If an invalid
         * sampling frequency is provided, an exception is thrown.
         *
         * @param[in] AdB Desired stopband attenuation for the halfband filter.
         * @param[in] fmax Maximum signal frequency present in the input.
         * @param[in] fs Sampling frequency (must be > 0).
         * @param[in] input Complex input samples.
         * @return std::vector<std::complex<double>> Interpolated output signal.
         *
         * @throws std::invalid_argument If `fs <= 0.0`.
         */
        std::vector<std::complex<double>> operator()(
            double AdB,
            double fmax,
            double fs,
            const std::vector<std::complex<double>>& input
        );

    private:
        /**
         * @brief Apply an FIR filter to a complex signal.
         *
         * Pads the input to avoid circular convolution artifacts, applies
         * the FIR filter using dsp::lfilter, and removes the padding.
         *
         * @param[in] b FIR filter coefficients.
         * @param[in] input Complex input samples.
         * @return std::vector<std::complex<double>> Filtered output.
         */
        std::vector<std::complex<double>> filter(
            const std::vector<double>& b,
            const std::vector<std::complex<double>>& input
        );

        /**
         * @brief Upsample a complex signal by a factor of two.
         *
         * Inserts a zero between each input sample, doubling the length.
         *
         * @param[in] input Complex input samples.
         * @return std::vector<std::complex<double>> Upsampled signal.
         */
        std::vector<std::complex<double>> upsample(
            const std::vector<std::complex<double>>& input
        );

        uint nSteps;
        HalfBand * halfband;
};
