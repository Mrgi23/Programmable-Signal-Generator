#ifndef DSP_H
#define DSP_H

#ifdef __cplusplus

#include <complex>
#include <vector>

namespace dsp {
    class ComplexFFT {
        private:
            std::vector<std::complex<double>> chirpZ(unsigned int N, const std::vector<std::complex<double>>& x);
            std::vector<std::complex<double>> dft(const std::vector<std::complex<double>>& x);
            std::vector<std::complex<double>> radix2(const std::vector<std::complex<double>>& x);
        public:
            ComplexFFT() {}
            ~ComplexFFT() {}

            template <typename T>
            std::vector<std::complex<double>> operator()(const std::vector<T>& x, unsigned int N = 0) {
                // Adjust the number of points for FFT, accordingly.
                unsigned int Nfft = N;
                if (!Nfft) { Nfft = x.size(); }

                // Cast input type to complex.
                std::vector<std::complex<double>> xfft;
                for (T sample : x) { xfft.push_back(static_cast<std::complex<double>>(sample)); }

                // Zero-padding, if necessary.
                if (xfft.size() < Nfft) { xfft.resize(Nfft, std::complex<double>(0.0, 0.0)); }

                // If number of FFT points is equal to the number of 2, compute FFT using radix2.
                if (!(Nfft & (Nfft - 1))) { return radix2(xfft); }

                // For smaller signals, compute FFT using DFT.
                if (x.size() <= Nfft && Nfft < 50) { return dft(xfft); }

                // Otherwise, compute FFT using chirpZ.
                return chirpZ(Nfft, xfft);
            }
    };
    std::vector<std::complex<double>> freqz(
        std::vector<double>& w,
        const std::vector<double>& b,
        unsigned int worN = 1024,
        double fs = 1.0
    );
    std::vector<std::complex<double>> lfilter(
        const std::vector<double>& b,
        const std::vector<std::complex<double>>& x
    );
    std::vector<double> remez(
        unsigned int numtaps,
        const std::vector<double>& bands,
        const std::vector<double>& desired,
        const std::vector<double>& weights = {},
        double fs = 1.0,
        unsigned int maxIter = 25,
        unsigned int gridDensity = 16
    );
}

#endif

#endif
