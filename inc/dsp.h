#ifndef DSP_H
#define DSP_H

#ifdef __cplusplus

#include <complex>
#include <type_traits>
#include <vector>

constexpr std::complex<double> complexJ = {0.0, 1.0};

template<typename T>
struct is_complex : std::false_type {};

template<typename T>
struct is_complex<std::complex<T>> : std::true_type {};

namespace dsp {
    class ComplexFFT {
        private:
            std::vector<std::complex<double>> chirpZ(uint N, const std::vector<std::complex<double>>& x);
            std::vector<std::complex<double>> dft(const std::vector<std::complex<double>>& x);
            std::vector<std::complex<double>> radix2(const std::vector<std::complex<double>>& x);
        public:
            ComplexFFT() {}
            ~ComplexFFT() {}

            template <typename T>
            std::vector<std::complex<double>> operator()(const std::vector<T>& x, uint N = 0) {
                // Adjust the number of points for FFT, accordingly.
                uint Nfft = N;
                if (!Nfft) { Nfft = x.size(); }

                // Cast input type to complex.
                std::vector<std::complex<double>> xfft;
                for (T sample : x) { xfft.push_back(static_cast<std::complex<double>>(sample)); }

                // Zero-padding, if necessary.
                if (xfft.size() < Nfft) { xfft.resize(Nfft, {0.0, 0.0}); }

                // If number of FFT points is equal to the number of 2, compute FFT using radix2.
                if (!(Nfft & (Nfft - 1))) { return radix2(xfft); }

                // For smaller signals, compute FFT using DFT.
                if (x.size() <= Nfft && Nfft < 50) { return dft(xfft); }

                // Otherwise, compute FFT using chirpZ.
                return chirpZ(Nfft, xfft);
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
    std::vector<double> remez(
        uint numtaps,
        const std::vector<double>& bands,
        const std::vector<double>& desired,
        const std::vector<double>& weights = {},
        double fs = 1.0,
        uint maxIter = 25,
        uint gridDensity = 16
    );
}

#endif

#endif
