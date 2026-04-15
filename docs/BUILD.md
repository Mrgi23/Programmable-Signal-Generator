# Programmable Signal Generator - Building Guide

## Overview
This project provides **two equivalent implementations**:
- **C++** – for high-performance native builds
- **Python** – for scripting and cross-platform workflows

Both implementations expose the same functionality, and you can build and run either independently.

## General Requirements (All Platforms)
## General Requirements
- **Git** – version control
- **Conan ≥ 2.1** – required for pulling pulling dependencies
- **CMake ≥ 4.1** – build configuration
- **GCC ≥ 12** – C++20 compatible compiler
- **Make** or **Ninja** – build system
- **Python ≥ 3.12** – required for Python implementation

## Install Dependencies by OS
o build and run the project, ensure you have the required tools installed:

### 1. Linux (Ubuntu/Debian)
```sh
sudo apt update && sudo apt install -y build-essential cmake gcc g++ git liblapack-dev libopenblas-dev make ninja-build pipx python3.12 python3.12-venv
sudo pipx ensure path

sudo pipx install conan
conan remote add conancenter https://center2.conan.io
conan remote add artifactory https://conan.mrgi23.com/artifactory/api/conan/Conan-Index
```

### 3. Windows
- Install [Git for Windows](https://git-scm.com/downloads)
- Install [Python3.12+](https://www.python.org/downloads/)
    - Ensure `python` and `pip` are added to the system `PATH`

## Build
### **Clone the Repository**
Clone the project using SSH:
```sh
git clone git@github.com:Mrgi23/Programmable-Signal-Generator.git
cd Programmable-Signal-Generator
```

### C++
Pull the dependencies using Conan, generate Makefiles with CMake, and compile:
```sh
conan install . --build=missing
cmake -G Ninja --preset conan-release
cmake --build --preset conan-release
```

### Python
Create and activate a virtual environment (optional), and install dependencies:

#### 1. Linux
```sh
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

#### 2. Windows
```sh
python -m venv venv
source venv\Scripts\activate
pip install -r requirements.txt
```

## Usage
The project provides two different entry points, depending on which implementation you want to run:

- **C++ executable** → `bin/signalgen`
- **Pure Python implementation** → `app/signalgen.py`

### C++ Executable
```sh
# From bin folder
./signalgen
```

### Python
```sh
# From root folder
PYTHONPATH=./src/python python app/signalgen.py
```

## Troubleshooting

| **Issue**          | **Possible Fix** |
|--------------------|------------------|
| `cmake` not found | Install using `sudo apt install cmake` |
| Compiler errors         | Ensure `g++` version is supports C++20 (`g++ --version`). |
| Python version mismatch | Run `python3.12` explicitly if needed. |

## Next Steps
For general project information, see the [README](../README.md).

For details on testing and validation methods used in this system, see [Testing & Validation](TESTING.md).
