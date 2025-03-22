import numpy as np
from signal_generator import SignalGenerator

if __name__ == "__main__":
    print("Digital signal from ./data/digitalSignal.txt, will be generated as analog.")

    while True:
        dac_mode = int(input("Please enter DAC reconstruction mode (0(NRZ) or 1(RF)): "))

        if dac_mode != 0 and dac_mode != 1:
            print("Invalid DAC reconstruction mode!")
            continue
        else:
            break

    # Signal parameters.
    data = np.loadtxt("./data/digitalSignal.txt", delimiter=",")
    digital = data[:, 0] + 1j * data[:, 1]
    N = len(digital)
    fs = 61.44e6
    f_max = 25e6

    # Interpolator parameters.
    n_steps = 4
    n_points = 8192
    A_dB = 60.0

    # CORDIC parameters.
    n_iter = 13
    f_res = 1.0
    if not dac_mode:
        f_shift = 368.64e6
    else:
        f_shift = 860.16e6

    # DAC parameters.
    if not dac_mode:
        mode = "NRZ"
    else:
        mode = "RF"
    n_nyquist = 4
    F_pass = 0.4
    error_dB = 0.025

    # Signal generator.
    gen = SignalGenerator(n_steps, n_points, n_iter, f_res)
    analog = gen(digital, fs, f_max, f_shift, mode, A_dB, n_nyquist, F_pass, error_dB)

    print("Writing...")

    # Write analog signal.
    with open("./data/analogSignal.txt", "w") as file:
        for sample in analog:
            file.write(f"{sample:.18f}\n")

    print("Analog signal has been writen in ./data/analogSignal.txt")
