#pragma once

#include "types.h"
#include "firFilter.h"
#include <string>

/**
 * @class DAC
 * @brief Digital‑to‑analog reconstruction using selectable output modes.
 *
 * This class implements a configurable digital‑to‑analog converter (DAC) model
 * that reconstructs a continuous‑time signal from discrete digital samples.
 * It supports multiple reconstruction modes (e.g., NRZ, RF) and applies
 * appropriate reconstruction kernels and inverse‑sinc compensation.
 *
 * The DAC internally uses an InverseSinc filter to correct for the inherent
 * sinc‑shaped frequency response of zero‑order hold reconstruction. The number
 * of points used for reconstruction is configurable.
 */
class DAC
{
    public:
        /**
         * @brief Construct a new DAC object.
         *
         * @param[in] nPoints Number of points used for inverse‑sinc correction.
         */
        DAC(uint nPoints = 8192);

        /**
         * @brief Destructor.
         *
         * Releases the internal InverseSinc filter.
         */
        ~DAC();

        /**
         * @brief Reconstruct an analog‑like waveform from digital samples.
         *
         * Initializes the DAC with a given number of reconstruction points and
         * allocates the internal InverseSinc filter used for frequency‑response
         * correction.
         *
         * @param[in] digital Input digital samples.
         * @param[in] mode Reconstruction mode ("NRZ" or "RF").
         * @param[in] nNyquist Number of Nyquist zones (must be even for RF mode).
         * @param[in] Fpass Passband edge for inverse‑sinc correction.
         * @param[in] errordB Allowed ripple/error in dB for inverse‑sinc design.
         * @return std::vector<double> Reconstructed analog‑like output samples.
         */
        std::vector<double> operator()(
            const std::vector<double>& digital,
            std::string mode,
            uint nNyquist = 4U,
            double Fpass = 0.4,
            double errordB = 0.025
        );

    private:
        /**
         * @brief Generate the reconstruction kernel for the selected mode.
         *
         * Applies the selected reconstruction mode (e.g., NRZ, RF), generates
         * the appropriate reconstruction kernel, applies inverse‑sinc
         * compensation, and returns the reconstructed output waveform.
         *
         * The RF mode requires an even number of Nyquist zones. If an invalid
         * mode is provided or the parameters are inconsistent, an exception is
         * thrown.
         *
         * @param[in] mode Reconstruction mode ("NRZ" or "RF").
         * @param[in] nNyquist Number of Nyquist zones.
         * @return std::vector<double> Reconstruction kernel.
         *
         * @throws std::invalid_argument If the mode is invalid.
         * @throws std::invalid_argument If RF mode is selected with an odd nNyquist.
         */
        std::vector<double> kernel(std::string mode, uint nNyquist);

        InverseSinc * inverseSinc;
};
