# Programmable Signal Generator

## Description

The **Programmable Signal Generator** is a fully software-simulated system designed for high-precision signal synthesis and modulation. It enables the generation of custom waveforms with fine frequency resolution, making it suitable for signal processing simulations, software-defined radio (SDR) testing, and general-purpose waveform generation.

### **Key Features**
- **Fully Software-Based** – No external hardware is required.
- **High-Precision Interpolation & Frequency Shifting** – Uses a CORDIC-based complex mixer for precise frequency shift.
- **Multi-Language Support** – Implemented in both **Python** and **C++** for flexible performance trade-offs.
- **Multiple Output Modes** – Supports NRZ and RF DAC simulation.
- **Ideal for Research & Education** – A valuable tool for algorithm development, DSP research, and software-based testing.

## Documentation
- [**System Architecture**](docs/ARCHITECTURE.md) – Structural breakdown of the system.
- [**Mathematical Background**](docs/MATH.md) – Theoretical and mathematical principles.
- [**Testing & Validation**](docs/TESTING.md) – Methods for verifying accuracy and performance.

## Setup
To install dependencies and run the project:

### **Clone the Repository**
Clone the project using SSH:
```bash
git clone git@gitlab.com:Mrgi23/programmable-signal-generator.git
cd programmable-signal-generator
```
### **Setup for Python**
Create a virtual environment (optional) and install dependencies:
```bash
python3 -m venv venv   # Optional: Create virtual environment
source venv/bin/activate  # Activate (Linux/macOS)
venv\Scripts\activate     # Activate (Windows)

pip install -r requirements.txt
```
### **Setup for C++**
Create a build/ directory, generate Makefiles with CMake, and compile:
```bash
mkdir -p build && cd build
cmake ..
make
```
### **Run the Program**
* **Python**
    ```bash
    python python/main.py
    ```
* **C++**
    ```bash
    ./build/signalGen.exe # Run (Linux/macOS)
    \build\signalGen.exe # Run (Windows)
    ```

## Authors and Acknowledgments
### Author
- [**Aleksa Markovic**](https://gitlab.com/Mrgi23) – Creator & Maintainer

### Acknowledgments
Special thanks to:
- [**School of Electrical Engineering, University of Belgrade**](https://www.etf.bg.ac.rs/) - For inspiring this project.
- The open-source community for valuable resources and support.

## License
This project is licensed under **CC0 1.0 Universal (Public Domain Dedication)**.

- **You are free to fork this repository and use it in your own project.**
- **You may NOT modify this repository directly.**
- **This project is in the public domain (no copyright), but any changes must be made in your own fork**.

For full legal terms, see the [LICENSE](LICENSE) file.

[![CC0](https://licensebuttons.net/p/zero/1.0/88x31.png)](https://creativecommons.org/publicdomain/zero/1.0/)

## Contributing
This repository does not accept direct contributions. However, you are free to fork and modify the code for your own use.

If you encounter any issues, feel free to open an issue, but please understand that fixes and feature requests may not be actively addressed.

## Next Steps
For a detailed breakdown of the system’s architecture, see the [System Architecture](docs/ARCHITECTURE.md).
