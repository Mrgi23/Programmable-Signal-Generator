#ifndef DSP_H
#define DSP_H

#ifdef __cplusplus

#include <algorithm>
#include <armadillo>
#include <complex>
#include <type_traits>
#include <vector>

template<typename T>
struct is_complex : std::false_type {};

template<typename T>
struct is_complex<std::complex<T>> : std::true_type {};

namespace dsp {
    class ComplexFFT {
        public:
            ComplexFFT() {}
            ~ComplexFFT() {}

            template <typename T>
            std::vector<std::complex<double>> fft(const std::vector<T>& x, uint N = 0) {
                // Adjust the number of points for FFT, accordingly.
                uint Nfft = N;
                if (!Nfft) { Nfft = x.size(); }

                // Compute FFT with Armadillo.
                arma::cx_vec armaXfft;
                if constexpr (std::is_same_v<T, double>) {
                    arma::vec armaX(x);
                    armaXfft = arma::fft(armaX, Nfft);
                }
                else if constexpr (std::is_same_v<T, std::complex<double>>) {
                    arma::cx_vec armaX(x);
                    armaXfft = arma::fft(armaX, Nfft);
                }

                // Cast Armadillo vector to standard complex vector.
                std::vector<std::complex<double>> xfft(armaXfft.begin(), armaXfft.end());
                return xfft;
            }

            template <typename T>
            std::vector<T> fftshift(const std::vector<T>& x) {
                // Define the shifted output.
                std::vector<T> shifted(x.size());

                // Calculate the pivot for the shift.
                uint middle = shifted.size() / 2;

                // Shift all elements around the pivot.
                for (uint i = 0; i < shifted.size(); i++) {
                    uint shift;
                    if (x.size() % 2) { shift = (i + middle + 1) % x.size(); }
                    else { shift = (i + middle) % x.size(); }
                    shifted[i] = x[shift];
                }
                return shifted;
            }
    };
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
        // Define the result type.
        using resultT = typename std::conditional<
            is_complex<T>::value,
            std::complex<double>,
            double
        >::type;

        // Calculate size of the kernel and the signal.
        uint nKernel = h.size();
        uint nSignal = x.size();

        // Calculate size of the convolution output.
        uint nConv = nSignal + nKernel - 1;

        // Define empty convolved signal.
        std::vector<resultT> y(nConv, static_cast<resultT>(0));

        // Perform the convolution between the input signal and the kernel.
        for (uint n = 0; n < nConv; n++) {
            for (uint k = 0; k < nKernel; k++) {
                // If k exceeds the current output index, no more contributions.
                if (k > n)
                    break;

                // Skip if (n - k) is outside the input signal range.
                if ((n - k) >= nSignal)
                    continue;

                y[n] += h[k] * x[n - k];
            }
        }
        return y;
    };
    std::vector<double> firls(
        uint numtaps,
        const std::vector<double>& bands,
        const std::vector<double>& desired,
        const std::vector<double>& weights = {},
        double fs = 2.0,
        uint gridSize = 1024
    );
    std::vector<std::complex<double>> freqz(
        std::vector<double>& w,
        const std::vector<double>& b,
        uint worN = 1024,
        double fs = 1.0
    );
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
        // Define the result type.
        using resultT = typename std::conditional<
            is_complex<T>::value,
            std::complex<double>,
            double
        >::type;

        // Calculate size of the filter and the signal.
        uint nFilter = b.size();
        uint nSignal = x.size();

        // Define empty filtered signal.
        std::vector<resultT> y(nSignal, static_cast<resultT>(0));

        // Filter the input signal with the filter coefficients.
        for (uint n = 0; n < nSignal; n++) {
            for (uint k = 0; k < nFilter && k <= n; k++) { y[n] += b[k] * x[n - k]; }
        }
        return y;
    };
}

#endif

#endif
