#ifndef UTILS_H
#define UTILS_H

#ifdef __cplusplus

#include <complex>
#include <string>
#include <vector>

namespace utils {
    template <typename T>
    std::vector<T> linspace(T start, T end, uint num, bool endpoint = true) {
        // Initialize linspace vector.
        std::vector<T> vec;

        // If number of points is 0, return empty vector.
        if (!num) { return vec; }

        // If number of points is 1, return start.
        if (num == 1) {
            vec.push_back(start);
            return vec;
        }

        // Fill in the linspace vector, with or without end element depending on endpoint.
        T step = endpoint ? (end - start) / static_cast<T>(num - 1) : (end - start) / static_cast<T>(num);
        for (uint i = 0; i < num; i++) { vec.push_back(start + i * step); }
        return vec; }; // linspace
    bool readFile(std::string path, std::vector<std::complex<double>>& signal);
    bool writeFile(std::string path, const std::vector<double>& signal);
}

#endif

#endif
