# Programmable Signal Generator - Testing & Validation Guide

## Overview
Testing strategies ensure the accuracy, reliability, and performance of the Programmable Signal Generator. Testing is divided into unit tests, integration tests, and coverage analysis, ensuring every component functions correctly in isolation and as part of the overall system.

## Testing Frameworks
### C++ Testing
- **Google Test (gtest), Google Mock (gmock)** – Testing framework for C++
- **LLVM `llvm-cov`** – Code coverage analysis

### Python Testing
- **pytest** – Testing framework for Python
- **pytest-cov** – Code coverage analysis

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

### Installation
#### 1. Linux/WSL
```sh
sudo apt update && sudo apt install -y git cmake clang lld make python3.10 python3.10-venv llvm lcov liblapack-dev libopenblas-dev

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

#### 2. macOS
```sh
brew install git cmake clang llvm make python@3.10 lcov

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Unit Testing
Unit testing ensures that individual components of the system function correctly in isolation. Each module is tested independently to verify expected behavior under different conditions.

### Scope of Unit Testing
- Validate **Interpolator**, **Complex Mixer**, **Sinc Compensation Filter**, and **Digital-to-Analog Converter**.
- Verify **data processing functions and signal transformations**.
- Verify **numerical stability, edge case handling, and computational performance**.

## Integration Testing
Integration tests verify that the various components of the Programmable Signal Generator work together correctly. These tests ensure that data flows properly between modules, that signal transformations are performed as expected, and that the overall system produces the correct output when given a known input.

### Scope of Integration Testing
- Verify core modules, **Interpolator**, **Complex Mixer**, **Sinc Compensation Filter**, and **Digital-to-Analog Converter interact correctly, ensuring data is passed seamlessly between them**.
- Ensure **complete processing chain produces the expected output**.
- Assess **performance and robustness of the overall system when all components are integrated**.

##  Running Tests
### C++ Tests
```sh
mkdir -p tests/build && cd tests/build
cmake ..
make
make unit # Unit Tests
make integration # Integration Tests
```

### Python Tests
```sh
source .venv/bin/activate
PYTHONPATH=./src/python pytest tests/python/unit # Unit Tests
PYTHONPATH=./src/python pytest tests/python/integration # Integration Tests
```

### Test File Location
| **Language**              | **Directory**              |
|---------------------------|----------------------------|
| C++ Unit Tests            | `tests/cpp/unit`           |
| C++ Integration Tests     | `tests/cpp/integration`    |
| Python Unit Tests         | `tests/python/unit`        |
| Python Integration Tests  | `tests/python/integration` |

For detailed test cases, see the corresponding test files in `tests/<language>/<type>`.

## **Code Coverage**
Code coverage ensures tests sufficiently exercise the codebase, identifying tested portions.

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
PYTHONPATH=./src/python pytest --cov=./ --cov-report=html:reports/htmlcov/python tests/python
```

### Coverage Report Locations
| **Language**       | **Directory**            |
|--------------------|--------------------------|
| C++ Unit Tests     | `reports/htmlcov/cpp`    |
| Python Unit Tests  | `reports/htmlcov/python` |

To view coverage report:
```sh
firefox <directory>/index.html
```

## CI/CD Integration
The testing framework is integrated with **GitLab CI/CD**, ensuring automated testing and coverage reporting on every merge request.

### CI/CD Pipeline Overview
The pipeline is structured into the following stages:
| **Stage**    | **Purpose** |
|--------------|-------------|
| **Setup**    | Builds and pushes a custom Docker image with all dependencies pre-installed (`cmake`, `clang`, `llvm`, `python3.10`, `googletest`, `googlemock`, `pytest`, etc.) |
| **Build**    | Compiles C++ tests |
| **Test**     | Runs both C++ (Google Test) and Python (pytest) tests, and generates coverage reports using `llvm-cov` and `pytest-cov` |
| **Deploy**   | Creates new tag based on the release version and publishes reports via **GitLab Pages** (only `main`) |

### When Does CI/CD Run?
The CI/CD pipeline is triggered in the following cases:
- On every commit and merge request to `develop` or `main` branches
- On manual pipeline execution from **GitLab UI**

**The pipeline is blocked if tests fail**, ensuring only validated code is merged.

### GitLab CI/CD Coverage Badge
![Coverage](https://gitlab.com/mrgi23/programmable-signal-generator/badges/main/coverage.svg)

### GitLab Pages (Published Reports)
Deployed via GitLab Pages, coverage reports are accessible at:
- [Programmable Signal Generator's Processing Pipeline](https://mrgi23.gitlab.io/programmable-signal-generator/index.html)
- [C++ Code Coverage Report](https://mrgi23.gitlab.io/programmable-signal-generator/cpp/index.html)
- [Python Code Coverage Report](https://mrgi23.gitlab.io/programmable-signal-generator/python/index.html)


## Next Steps
For instructions on building and running the project, see the [Building Guide](BUILD.md).

For a more in-depth look into the mathematical principles behind each component, see the [Mathematical Background](MATH.md).
