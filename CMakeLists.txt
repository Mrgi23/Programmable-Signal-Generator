cmake_minimum_required(VERSION 3.10)

# Declare the project.
project(ProgrammableSignalGenerator)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find dependencies
find_package(armadillo REQUIRED)
find_package(GTest REQUIRED)
find_package(liquid-dsp REQUIRED)

# Set ROOT_DIR
get_filename_component(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}" ABSOLUTE)

# Include project headers
include_directories("${ROOT_DIR}/inc")

# Collect all source files from src/ directory
file(GLOB CXX_SOURCES "${ROOT_DIR}/src/cpp/*.cpp")

# Create a single test executable
if(WIN32)
  add_executable(signalGen.exe ${CXX_SOURCES})
else()
  add_executable(signalGen ${CXX_SOURCES})
endif()

# Link with Armadillo and LiquidDSP
target_link_libraries(signalGen armadillo::armadillo liquid-dsp::liquid-dsp)
