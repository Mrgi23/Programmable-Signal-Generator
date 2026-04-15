#pragma once

#include "types.h"
#include <algorithm>
#include <armadillo>
#include <complex>
#include <type_traits>
#include <vector>

/**
 * @brief Type trait for detecting std::complex<T>.
 *
 * Evaluates to true if T is a specialization of std::complex, otherwise false.
 */
template<typename T>
struct is_complex : std::false_type {};

template<typename T>
struct is_complex<std::complex<T>> : std::true_type {};

/**
 * @namespace dsp
 * @brief Core digital signal processing utilities and algorithms.
 *
 * This namespace contains foundational DSP components used throughout the
 * signal‑generation pipeline. It provides:
 *
 * - FFT utilities (`ComplexFFT`): Real and complex FFT computation and
 *   frequency‑domain reordering (fftshift).
 *
 * - Convolution and filtering (`convolve`, `lfilter`): Generic FIR
 *   convolution and direct‑form filtering for both real and complex signals.
 *
 * - Filter design algorithms (`firls`, `freqz`): Least‑squares FIR filter
 *   design and frequency‑response evaluation.
 *
 * The functions and classes in this namespace are designed to be lightweight,
 * modular, and compatible with both real and complex data types. They form the
 * mathematical backbone for higher‑level DSP components such as interpolators,
 * mixers, and DAC reconstruction filters.
 */
namespace dsp
{
    /**
     * @class ComplexFFT
     * @brief FFT utilities for real and complex input vectors.
     *
     * Provides a unified interface for computing FFTs on real-valued or
     * complex-valued vectors using Armadillo’s FFT backend. Also includes
     * an FFT-shift operation that reorders frequency bins into centered form.
     */
    class ComplexFFT
    {
        public:
            /**
             * @brief Construct a new ComplexFFT object.
             *
             */
            ComplexFFT();

            /**
             * @brief Destructor.
             *
             * Default cleanup.
             */
            ~ComplexFFT();

            /**
             * @brief Compute the FFT of a real or complex input vector.
             *
             * Compute an N‑point FFT. If N is zero, the FFT length defaults to the input vector size.
             *
             * Supported input types:
             * - double
             * - std::complex<double>
             *
             * @tparam T Input sample type.
             * @param[in] x Input vector.
             * @param[in] N Optional FFT length (0 → use x.size()).
             * @return std::vector<std::complex<double>> FFT output.
             */
            template <typename T>
            std::vector<std::complex<double>> fft(const std::vector<T>& x, uint N = 0)
            {
                uint Nfft = N;
                if (!Nfft)
                    Nfft = x.size();

                arma::cx_vec armaXfft;
                if constexpr (std::is_same_v<T, double>)
                {
                    arma::vec armaX(x);
                    armaXfft = arma::fft(armaX, Nfft);
                }
                else if constexpr (std::is_same_v<T, std::complex<double>>)
                {
                    arma::cx_vec armaX(x);
                    armaXfft = arma::fft(armaX, Nfft);
                }

                std::vector<std::complex<double>> xfft(armaXfft.begin(), armaXfft.end());
                return xfft;
            }

            /**
             * @brief Perform FFT-shift on a vector.
             *
             * Reorders the input vector so that the zero-frequency component
             * is centered. For even-length vectors, the shift is symmetric.
             * For odd-length vectors, the center index is offset by one.
             *
             * @tparam T Element type.
             * @param[in] x Input vector.
             * @return std::vector<T> Shifted output vector.
             */
            template <typename T>
            std::vector<T> fftshift(const std::vector<T>& x)
            {
                std::vector<T> shifted(x.size());
                uint middle = shifted.size() / 2;
                for (uint i = 0; i < shifted.size(); i++)
                {
                    uint shift;
                    if (x.size() % 2) { shift = (i + middle + 1) % x.size(); }
                    else { shift = (i + middle) % x.size(); }
                    shifted[i] = x[shift];
                }
                return shifted;
            }
    };

    /**
     * @brief Convolution of a real kernel with real or complex input.
     *
     * Computes linear convolution between a real-valued FIR kernel `h`
     * and a signal `x`, which may be real or complex. The output type
     * automatically matches the input type.
     *
     * @tparam T Input sample type (double or std::complex<double>).
     * @param[in] h FIR kernel coefficients.
     * @param[in] x Input signal.
     * @return Convolution result with appropriate type.
     */
    template <typename T>
    auto convolve(
        const std::vector<double>& h,
        const std::vector<T>& x
    ) -> typename std::conditional<
            is_complex<T>::value,
            std::vector<std::complex<double>>,
            std::vector<double>
        >::type
    {
        using resultT = typename std::conditional<
            is_complex<T>::value,
            std::complex<double>,
            double
        >::type;

        uint nKernel = h.size();
        uint nSignal = x.size();
        uint nConv = nSignal + nKernel - 1;

        std::vector<resultT> y(nConv, static_cast<resultT>(0));
        for (uint n = 0; n < nConv; n++)
            for (uint k = 0; k < nKernel; k++)
            {
                if (k > n)
                    break;

                if ((n - k) >= nSignal)
                    continue;

                y[n] += h[k] * x[n - k];
            }
        return y;
    };

    /**
     * @brief Least-squares FIR filter design.
     *
     * Designs a linear-phase FIR filter using the least-squares method.
     * Band edges and desired responses are specified in pairs. Optional
     * weighting allows emphasizing certain bands.
     *
     * @param[in] numtaps Number of filter taps (must be odd).
     * @param[in] bands Frequency band edges (even-length vector).
     * @param[in] desired Desired response at each band edge.
     * @param[in] weights Optional per-band weights.
     * @param[in] fs Sampling frequency.
     * @param[in] gridSize Frequency grid resolution.
     * @return std::vector<double> FIR filter coefficients.
     *
     * @throws std::invalid_argument If numtaps is even.
     * @throws std::invalid_argument If bands has odd length.
     * @throws std::invalid_argument If desired size mismatches bands.
     * @throws std::invalid_argument If weights size mismatches band count.
     * @throws std::invalid_argument If fs <= 0.
     * @throws std::invalid_argument If band edges fall outside [0, 1].
     */
    std::vector<double> firls(
        uint numtaps,
        const std::vector<double>& bands,
        const std::vector<double>& desired,
        const std::vector<double>& weights = {},
        double fs = 2.0,
        uint gridSize = 1024
    );

    /**
     * @brief Compute the frequency response of an FIR filter.
     *
     * Evaluates the complex frequency response of FIR coefficients `b`
     * at `worN` frequency points. Frequencies are returned in `w`.
     *
     * @param[in,out] w Frequency vector (resized if needed).
     * @param[in] b FIR filter coefficients.
     * @param[in] worN Number of frequency points.
     * @param[in] fs Sampling frequency.
     * @return std::vector<std::complex<double>> Frequency response.
     */
    std::vector<std::complex<double>> freqz(
        std::vector<double>& w,
        const std::vector<double>& b,
        uint worN = 1024,
        double fs = 1.0
    );

    /**
     * @brief Apply an FIR filter to a real or complex signal.
     *
     * Computes the direct-form FIR filtering operation:
     * y[n] = Σ b[k] * x[n − k]
     *
     * The output type matches the input type.
     *
     * @tparam T Input sample type.
     * @param[in] b FIR filter coefficients.
     * @param[in] x Input signal.
     * @return Filtered output signal.
     */
    template <typename T>
    auto lfilter(
        const std::vector<double>& b,
        const std::vector<T>& x
    ) -> typename std::conditional<
            is_complex<T>::value,
            std::vector<std::complex<double>>,
            std::vector<double>
        >::type
    {
        using resultT = typename std::conditional<
            is_complex<T>::value,
            std::complex<double>,
            double
        >::type;

        uint nFilter = b.size();
        uint nSignal = x.size();

        std::vector<resultT> y(nSignal, static_cast<resultT>(0));
        for (uint n = 0; n < nSignal; n++)
            for (uint k = 0; k < nFilter && k <= n; k++)
                y[n] += b[k] * x[n - k];
        return y;
    };
}
