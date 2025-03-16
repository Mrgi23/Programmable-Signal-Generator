#include <iostream>
#include "utils.h"
#include "signalGenerator.h"

using namespace std;

int main() {
    std::cout << "Digital signal from ./data/digitalSignal.txt, will be generated as analog." << std::endl;

    int dacMode = -1;
    while (true) {
        cout << "Please enter DAC reconstruction mode (0(NRZ) or 1(RF)): ";
        cin >> dacMode;

        if (dacMode != 0 && dacMode != 1) {
            cout << "Invalid DAC reconstruction mode!" << endl;
            continue;
        }
        else { break; }
    }

    // Signal parameters.
    string inPath = "../data/digitalSignal.txt";
    vector<complex<double>> digital;
    utils::readFile(inPath, digital);
    uint N = digital.size();
    double fs = 61.44e6;
    double fmax = 25e6;

    // Interpolator parameters.
    uint nSteps = 4;
    uint nPoints = 8192;
    double AdB = 60.0;

    // CORDIC parameters.
    uint nIter = 13;
    double fres = 1.0;
    double fshift;
    if (!dacMode) { fshift = 368.64e6; }
    else { fshift = 860.16e6; }

    // DAC parameters.
    string mode;
    if (!dacMode) { mode = "NRZ"; }
    else { mode = "RF"; }
    uint nNyquist = 4;
    double Fpass = 0.4;
    double errordB = 0.025;

    // Signal generator.
    SignalGenerator gen(nSteps, nPoints, nIter, fres);
    vector<double> analog = gen(digital, fs, fmax, fshift, mode, AdB, nNyquist, Fpass, errordB);

    cout << "Writing..." << endl;

    // Write analog signal.
    string outPath = "../data/analogSignal.txt";
    utils::writeFile(outPath, analog);

    cout << "Analog signal has been writen in ./data/analogSignal.txt" << endl;

    return 0;
}