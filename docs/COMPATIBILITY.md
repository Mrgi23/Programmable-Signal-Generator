# Programmable Signal Generator - Testing Framework Compatibility

## Compatibility Issues
- The **Programmable Signal Generator** is designed to run tests on **Linux** and **macOS**.
- **On Windows,tests must be run inside Windows Subsystem for Linux (WSL).**

### Fully Supported Platforms
- **Linux (Ubuntu/Debian)**
- **macOS (Intel & Apple Silicon)**
- **Windows (via WSL)**

### Not Supported
- **Native Windows (MSVC, MinGW)**
  - `llvm-cov` does not work natively on Windows.
  - Google Test (`gtest`) requires manual setup on naitve Windows.
  - Profiling tools like (`Valgrind`) are not available on Windows.

## Installing WSL with Ubuntu
Open PowerShell as Administrator and run
```powershell
wsl --install -d Ubuntu
```
For Windows 10 users, ensure you have **WSL 2** installed:
```powershell
wsl --set-default-version 2
```

## Next Steps
For details on testing and validation methods used in this system, see [Testing & Validation](TESTING.md).
