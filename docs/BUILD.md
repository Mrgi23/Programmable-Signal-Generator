# Programmable Signal Generator - Building Guide

## General Requirements (All Platforms)
- **Git** - For cloning the repository
- **CMake (≥ 3.10)** - For configuring the C++ build
- **C++ Compiler supporting C++20**
    - **Linux:** `g++ (≥ 10)` or `clang++ (≥ 11)`
    - **macOS:** `clang++` (via Xcode command-line tools) or install `g++` via Homebrew
    - **Windows:** `MSVC (Visual Studio 2019+)` or `MinGW-w64 (≥ 10)` (for GCC)
- **Make or Ninja** - For compiling the C++ code
- **Python 3.10+** - Required for the Python part of the project

## Install Dependencies by OS
To install dependencies and run the project, ensure you have the required tools installed:

### 1. Linux (Ubuntu/Debian)
```sh
sudo apt update && sudo apt install -y git cmake g++ make python3.10 python3.10-venv
```
For `clang` instead of `gcc`:
```sh
sudo apt install -y clang lld
```

### 2. macOS (via Homebrew)
```sh
brew install git cmake make gcc python@3.10
```
To use `GNU Make` instead of `BSD Make`:
```sh
alias make=gmake
brew install make
```

### 3. Windows
- Install [Git for Windows](https://git-scm.com/downloads)
- Install [CMake](https://cmake.org/download/)
- Choose a C++ Compiler:
    - `Microsoft Visual Studio (MSVC)`
        - Install Visual Studio 2019 or later with the C++ CMake tools
        - Open "x64 Native Tools Command Prompt for VS" before running commands
    - `MinGW-w64 (GCC)`
        - Install via MSYS2:
            ```sh
            pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-make
            ```
        - Add `C:\msys64\mingw64\bin` to your system `PATH`
- Install [Python3.10+](https://www.python.org/downloads/)
    - Ensure `python` and `pip` are added to the system `PATH`

## Build
### **Clone the Repository**
Clone the project using SSH:
```sh
git clone git@gitlab.com:Mrgi23/programmable-signal-generator.git
cd programmable-signal-generator
```

### Setup for C++
Create a build/ directory, generate Makefiles with CMake, and compile:

#### 1. Linux/macOS
```sh
mkdir -p build && cd build
cmake ..
make
```

#### 2. Windows
- `Microsoft Visual Studio (MSVC)`:
    ```sh
    mkdir build && cd build
    cmake -G "Visual Studio 17 2022" ..
    cmake --build . --config Release
    ```
- `MinGW-w64 (GCC)`:
    ```sh
    mkdir build && cd build
    cmake -G "MinGW Makefiles" ..
    mingw32-make
    ```

### Setup for Python
Create and activate a virtual environment (optional), and install dependencies:

#### 1. Linux/macOS
```sh
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

#### 2. Windows
```sh
python3 -m venv venv
source venv\Scripts\activate
pip install -r requirements.txt
```

## Run the Program
### C++
#### 1. Linux/macOS
```sh
./build/signalGen.exe
```

#### 2. Windows
- `Microsoft Visual Studio (MSVC)`:
    ```sh
    build\Release\signalGen.exe
    ```
- `MinGW-w64 (GCC)`:
    ```sh
    build\signalGen.exe
    ```

### Python
#### 1. Linux/macOS
```sh
python src/python/main.py
```

#### 2. Windows
```sh
python src\python\main.py
```

## Troubleshooting

| **Issue**          | **Possible Fix** |
|--------------------|------------------|
| `cmake` not found | Install using `sudo apt install cmake` (Linux) or `brew install cmake` (macOS), or download from [cmake.org](https://cmake.org/download/) (Windows). |
| `make` not found | Use `gmake` on macOS, `mingw32-make` on Windows (MinGW), or `cmake --build .` for MSVC. |
| Compiler errors | Ensure your compiler supports C++20 (`g++ --version`, `clang++ --version`, or `cl.exe`). |
| Python version mismatch | Run `python3.10` explicitly if needed. |

## Next Steps
For general project information, see the [README](../README.md).

For details on testing and validation methods used in this system, see [Testing & Validation](TESTING.md).
