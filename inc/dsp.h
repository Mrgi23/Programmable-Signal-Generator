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

            std::vector<std::complex<double>> operator()(const std::vector<std::complex<double>>& x, unsigned int N = 0);
    };
    std::vector<std::complex<double>> freqz(
        std::vector<double>& w,
        const std::vector<double>& b,
        unsigned int worN = 1024,
        double fs = 1
    );
    std::vector<std::complex<double>> lfilter(
        const std::vector<double>& b,
        const std::vector<std::complex<double>>& x
    );
    std::vector<double> remez(
        unsigned int numtaps,
        const std::vector<double>& bands,
        const std::vector<double>& desired,
        const std::vector<double>& weight,
        double fs = 1.0,
        unsigned int maxIter = 25,
        unsigned int gridDensity = 16
    );
}

#endif

#endif
