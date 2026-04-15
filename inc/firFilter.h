#pragma once

#include "types.h"
#include <vector>

/**
 * @class FIR
 * @brief Base class for FIR‑based reconstruction and correction filters.
 *
 * Provides a common interface and storage for the number of filter points
 * used by derived FIR‑design classes. This class is not intended to be used
 * directly; it serves as a foundation for specialized FIR designs such as
 * inverse‑sinc correction and halfband filters.
 */
class FIR
{
    public:
        /**
         * @brief Construct a new FIR base object.
         *
         * @param[in] nPoints Number of frequency‑domain evaluation points.
         */
        FIR(uint nPoints = 8192);

        /**
         * @brief Virtual destructor.
         *
         * Allows safe destruction through base‑class pointers.
         */
        virtual ~FIR();

    protected:
        int nPoints;
};

/**
 * @class InverseSinc
 * @brief Designs an inverse‑sinc correction filter for DAC reconstruction.
 *
 * This class computes an FIR filter that compensates for the sinc‑shaped
 * frequency response of zero‑order hold DAC reconstruction. The filter is
 * designed iteratively using least‑squares FIR design until the passband
 * error meets the specified tolerance.
 */
class InverseSinc : public FIR
{
    public:
        /**
         * @brief Construct a new InverseSinc filter designer.
         *
         * @param[in] nPoints Number of frequency‑domain evaluation points.
         */
        InverseSinc(uint nPoints = 8192);

        /**
         * @brief Destructor.
         *
         * Default cleanup.
         */
        ~InverseSinc() override;

        /**
         * @brief Design an inverse‑sinc FIR filter.
         *
         * Generates a filter that approximates the inverse of the sinc response
         * over the passband `[0, Fpass]`. The design iteratively increases the
         * filter order until the magnitude response meets the allowed ripple
         * specified by `errordB`.
         *
         * @param[in] Fpass Passband edge (0.0 < Fpass < 0.5).
         * @param[in] errordB Allowed passband ripple in dB.
         * @param[in] nSpec Number of specification points for the target response.
         * @return std::vector<double> FIR filter coefficients.
         *
         * @throws std::invalid_argument If Fpass is outside (0, 0.5).
         * @throws std::length_error If the required filter order exceeds limits.
         */
        virtual std::vector<double> operator()(double Fpass, double errordB, uint nSpec = 16);
};

/**
 * @class HalfBand
 * @brief Designs a halfband low‑pass FIR filter.
 *
 * This class generates a halfband filter with a narrow transition band and
 * symmetric coefficients. The filter is designed using the Parks–McClellan
 * algorithm (via liquid‑dsp) and validated against the desired passband
 * ripple. If the design does not meet specifications, the filter order is
 * increased iteratively.
 */
class HalfBand : public FIR
{
    public:
        /**
         * @brief Construct a new HalfBand filter designer.
         *
         * @param[in] nPoints Number of frequency‑domain evaluation points.
         */
        HalfBand(uint nPoints = 8192);

        /**
         * @brief Destructor.
         *
         * Default cleanup.
         */
        ~HalfBand() override;

        /**
         * @brief Design a halfband FIR filter.
         *
         * Creates a halfband low‑pass filter with passband edge `Fpass` and
         * attenuation `AdB`. The design uses liquid‑dsp’s Parks–McClellan
         * implementation and validates the resulting frequency response.
         *
         * @param[in] AdB Desired stopband attenuation in dB.
         * @param[in] Fpass Passband edge (0.0 < Fpass < 0.25).
         * @return std::vector<double> Halfband FIR coefficients.
         *
         * @throws std::invalid_argument If Fpass is outside (0, 0.25).
         * @throws std::length_error If the required filter order exceeds limits.
         */
        virtual std::vector<double> operator()(double AdB, double Fpass);
};
