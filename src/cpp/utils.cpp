#include "utils.h"
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace std;

namespace utils
{
    bool readFile(string path, vector<complex<double>>& signal)
    {
        ifstream file(path);
        if (!file.is_open())
            return false;

        signal.clear();

        string line;
        while (getline(file, line)) {
            istringstream iss(line);
            double real, imag;
            char comma;

            iss >> real;
            iss >> comma;
            iss >> imag;
            signal.push_back({real, imag});
        }

        file.close();
        return true;
    }

    bool writeFile(string path, const vector<double>& signal)
    {
        ofstream file(path);

        ostringstream oss;
        oss << fixed << setprecision(18);
        for (double sample : signal)
            oss << sample << "\n";
        file << oss.str();

        file.close();
        return true;
    }
}
