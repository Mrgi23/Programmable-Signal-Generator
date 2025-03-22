#include <fstream>
#include <iomanip>
#include <sstream>
#include "utils.h"

using namespace std;

namespace utils {
    bool readFile(string path, vector<complex<double>>& signal) {
        // Open file.
        ifstream file(path);
        if (!file.is_open()) { return false; }

        // Clear the vector if it is not empty.
        signal.clear();

        // Read each line.
        string line;
        while (getline(file, line)) {
            istringstream iss(line);
            double real, imag;
            char comma;

            // Extract the real and the imaginary part of the number.
            iss >> real;
            iss >> comma;
            iss >> imag;

            // Fill in the signal.
            signal.push_back({real, imag});
        }
        // Close file.
        file.close();
        return true;
    }

    bool writeFile(string path, const vector<double>& signal) {
        // Open file.
        ofstream file(path);

        // Write each sample in the new line.
        ostringstream oss;
        oss << fixed << setprecision(18);
        for (double sample : signal) { oss << sample << "\n"; }
        file << oss.str();

        // Close the file.
        file.close();
        return true;
    }
}
