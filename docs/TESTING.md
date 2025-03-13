# Programmable Signal Generator - Testing & Validation Guide

## Overview
Testing strategies ensure the accuracy, reliability, and performance of the Programmable Signal Generator. Testing is divided into unit tests, integration tests, and coverage analysis, ensuring every component functions correctly in isolation and as part of the overall system.

## Testing Frameworks
### C++ Testing
- **Google Test (gtest), Google Mock (gmock)** – Testing framework for C++
- **LLVM `llvm-cov`** – Code coverage analysis
- **Valgrind** (optional) – Memory leak detection (Linux/macOS only)

### Python Testing
- **pytest** – Testing framework for Python
- **pytest-cov** – Code coverage analysis
- **mypy** – Static type checking for Python

### Compatibility Notice
For compatibility details, see the [Testing Framework Compatibility](COMPATIBILITY.md)

## Installing Dependencies
Ensure that all required dependencies are installed.

### C++ Dependencies
- **CMake** – Build system for compiling tests
- **LLVM (`llvm-cov`)** – Required for C++ code coverage (**Only works with Clang**)
- **Make or Ninja** – For compiling tests

**C++ code coverage is only supported with Clang**
- If using **GCC or MSVC**, code coverage will not be available.
- Windows users must run tests inside **WSL with Clang/LLVM**.

### Python Dependencies
- **pytest** – Python unit testing framework
- **pytest-cov** –  Python code coverage tool
- **mypy** – Static type checking

### Installation
#### 1. Linux/WSL
```sh
sudo apt update && sudo apt install -y git cmake clang lld make python3.10 python3.10-venv llvm lcov

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

#### 2. macOS
```sh
brew install git cmake clang llvm make python@3.10 lcov

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Unit Testing
Unit testing ensures that individual components of the system function correctly in isolation. Each module is tested independently to verify expected behavior under different conditions.

### Scope of Unit Testing
- Validate **Interpolator**, **Complex Mixer**, **Sinc Compensation Filter**, and **Digital-to-Analog Converter**.
- Verify **data processing functions and signal transformations**.
- Verify **numerical stability, edge case handling, and computational performance**.

###  Running Unit Tests
#### C++ Tests
```sh
mkdir -p tests/build && cd tests/build
cmake ..
make
./unitTests
```

#### Python Tests
```sh
source .venv/bin/activate
PYTHONPATH=./src/python pytest tests/python/unit
```

### Test File Location
| **Language**       | **Directory**    |
|--------------------|------------------|
| C++ Unit Tests     | `tests/cpp/unit`      |
| Python Unit Tests  | `tests/python/unit`   |

For detailed test cases, see the corresponding test files in `tests/`.

## **Code Coverage**
Code coverage ensures tests sufficiently exercise the codebase, identifying untested portions.

### Generating Coverage Reports
#### C++ Code Coverage
C++ code coverage is generated using **LLVM's `llvm-cov`**, which only works with **Clang**.
**MSVC and MinGW are not supported for code coverage.**
- Linux/macOS: **Native support with Clang**
- Windows: **Must use WSL with Clang/LLVM**
```sh
mkdir -p tests/build && cd tests/build
cmake ..
make
make coverage
```
#### Python Code Coverage
```sh
source .venv/bin/activate
PYTHONPATH=./src/python pytest --cov=./ --cov-report=html:htmlcov/python tests/python
```

### Coverage Report Locations
| **Language**       | **Directory**    |
|--------------------|------------------|
| C++ Unit Tests     | `htmlcov/cpp`      |
| Python Unit Tests  | `htmlcov/python`   |

To view coverage report:
```sh
firefox directory/index.html
```

## CI/CD Integration
The testing framework is integrated with **GitLab CI/CD**, ensuring automated testing and coverage reporting on every merge request.

### CI/CD Pipeline Overview
The pipeline is structured into the following stages:
| **Stage**    | **Purpose** |
|--------------|-------------|
| **Setup**    | Installs dependencies (CMake, Clang, Python, Google Test, pytest) |
| **Build**    | Compiles C++ tests |
| **Test**     | Runs C++ (Google Test) and Python (pytest) tests |
| **Coverage** | Generates test coverage reports for C++ (`llvm-cov`) and Python (`pytest-cov`) |
| **Deploy**   | Publishes coverage reports via **GitLab Pages** |

### When Does CI/CD Run?
The CI/CD pipeline is triggered in the following cases:
- On every merge request to `develop` or `main` branches
- On manual pipeline execution from **GitLab UI**

**The pipeline is blocked if tests fail**, ensuring only validated code is merged.

### GitLab CI/CD Coverage Badge
![Coverage](https://gitlab.com/mrgi23/programmable-signal-generator/badges/main/coverage.svg)

### GitLab Pages (Published Reports)
Deployed via GitLab Pages, coverage reports are accessible at:
- [C++ Code Coverage Report](https://mrgi23.gitlab.io/programmable-signal-generator/cpp/index.html)
- [Python Code Coverage Report](https://mrgi23.gitlab.io/programmable-signal-generator/python/index.html)


## Next Steps
For instructions on building and running the project, see the [Building Guide](BUILD.md).

For a more in-depth look into the mathematical principles behind each component, see the [Mathematical Background](MATH.md).
