#pragma once

#include "types.h"
#include "interpolator.h"
#include "complexMixer.h"
#include "dac.h"
#include <string>

/**
 * @class SignalGenerator
 * @brief Full digital signal generation chain: interpolation → frequency shift → DAC reconstruction.
 *
 * This class implements a complete multirate digital signal generation pipeline.
 * It takes a complex baseband signal, increases its sampling rate through multistage
 * interpolation, applies a frequency shift using a complex mixer, and finally performs
 * digital‑to‑analog reconstruction using a selectable DAC mode.
 *
 * The generator combines three major DSP components:
 *
 * - Interpolator: Multistage ×2 upsampling using halfband filters.
 *
 * - ComplexMixer: Frequency translation using CORDIC/NCO‑based mixing.
 *
 * - DAC: Reconstruction filtering with selectable NRZ/RF modes and inverse‑sinc correction.
 *
 * This class provides a single high‑level operator() that executes the entire chain.
 */
class SignalGenerator {
    public:
        /**
         * @brief Construct a new SignalGenerator.
         *
         * @param[in] nSteps Number of ×2 interpolation stages.
         * @param[in] nPoints Number of frequency‑domain points for filter design.
         * @param[in] nIter Number of CORDIC iterations for the mixer.
         * @param[in] fres Frequency resolution for the mixer’s NCO.
         */
        SignalGenerator(uint nSteps = 4, uint nPoints = 8192, uint nIter = 13, double fres = 1.0);

        /**
         * @brief Destructor.
         *
         * Frees all internal DSP components.
         */
        ~SignalGenerator();

        /**
         * @brief Generate an analog‑like waveform from a complex baseband signal.
         *
         * The processing chain is:
         *
         * 1. Interpolation: The input signal is upsampled by `2^nSteps` using halfband filters.
         *
         * 2. Frequency shift: The interpolated I/Q components are mixed with a complex exponential
         *    at frequency `fshift`, using the ComplexMixer.
         *
         * 3. DAC reconstruction: The shifted signal is passed through the DAC model, which applies
         *    the selected reconstruction mode (`mode`) and inverse‑sinc correction.
         *
         * The sampling frequency `fs` must be strictly positive. If invalid, an exception is thrown.
         *
         * @param[in] signal Complex baseband input samples.
         * @param[in] fs Original sampling frequency (must be > 0).
         * @param[in] fmax Maximum signal frequency present in the input.
         * @param[in] fshift Frequency shift to apply after interpolation.
         * @param[in] mode DAC reconstruction mode ("NRZ" or "RF").
         * @param[in] AdB Stopband attenuation for halfband filters.
         * @param[in] nNyquist Number of Nyquist zones for DAC reconstruction.
         * @param[in] Fpass Passband edge for inverse‑sinc correction.
         * @param[in] errordB Allowed ripple for inverse‑sinc correction.
         * @return std::vector<double> Final reconstructed analog‑like waveform.
         */
        std::vector<double> operator()(
            const std::vector<std::complex<double>>& signal,
            double fs,
            double fmax,
            double fshift,
            std::string mode,
            double AdB = 60.0,
            uint nNyquist = 4U,
            double Fpass = 0.4,
            double errordB = 0.025
        );

    private:
        Interpolator * interpolator;
        ComplexMixer * complexMixer;
        DAC * dac;
};
