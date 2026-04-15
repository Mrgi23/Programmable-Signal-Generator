#pragma once

#include "types.h"
#include <complex>
#include <string>
#include <vector>

/**
 * @namespace utils
 * @brief Utility helpers for numerical operations and file I/O.
 *
 * This namespace provides small, reusable helper functions used throughout
 * the DSP pipeline. It includes:
 *
 * - A templated `linspace` generator for producing evenly spaced values.
 *
 * - File readers/writers for complex and real-valued signals.
 *
 * These utilities are lightweight and do not depend on external DSP libraries.
 */
namespace utils
{
    /**
     * @brief Generate a linearly spaced vector.
     *
     * Produces `num` evenly spaced values between `start` and `end`. If
     * `endpoint` is true, the last value equals `end`; otherwise, the range
     * is divided into `num` equal intervals without including the endpoint.
     *
     * @tparam T Numeric type (e.g., double, float).
     * @param[in] start Starting value.
     * @param[in] end Ending value.
     * @param[in] num Number of points to generate.
     * @param[in] endpoint Whether to include the endpoint.
     * @return std::vector<T> Linearly spaced values.
     */
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

        T step = endpoint ? (end - start) / static_cast<T>(num - 1)
                          : (end - start) / static_cast<T>(num);

        for (uint i = 0; i < num; i++)
            vec.push_back(start + i * step);

        return vec; };

    /**
     * @brief Read a complex-valued signal from a text file.
     *
     * The file must contain one complex sample per line in the format real,imag
     *
     * Whitespace is ignored. Each parsed pair is appended to the output vector.
     *
     * @param[in] path Path to the input file.
     * @param[out] signal Parsed complex-valued samples.
     * @return true If the file was successfully opened and parsed.
     * @return false If the file could not be opened.
     */
    bool readFile(std::string path, std::vector<std::complex<double>>& signal);

    /**
     * @brief Write a real-valued signal to a text file.
     *
     * Writes one sample per line with fixed high precision (18 digits).
     *
     * @param[in] path Output file path.
     * @param[in] signal Real-valued samples to write.
     * @return true Always returns true if the file could be opened.
     */
    bool writeFile(std::string path, const std::vector<double>& signal);
}
