#pragma once

#include "types.h"
#include <complex>
#include <string>
#include <vector>

namespace utils
{
    template <typename T>
    std::vector<T> linspace(T start, T end, uint num, bool endpoint = true)
    {
        std::vector<T> vec;
        if (!num)
            return vec;
        else if (num == 1)
        {
            vec.push_back(start);
            return vec;
        }

        T step = endpoint ? (end - start) / static_cast<T>(num - 1) : (end - start) / static_cast<T>(num);
        for (uint i = 0; i < num; i++) { vec.push_back(start + i * step); }
        return vec; };

    bool readFile(std::string path, std::vector<std::complex<double>>& signal);
    bool writeFile(std::string path, const std::vector<double>& signal);
}
