cmake_minimum_required(VERSION 3.10)

# Declare the project.
project(ProgrammableSignalGenerator)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Includes
include(FetchContent)

# Fetch Armadillo
FetchContent_Declare(
  armadillo
  GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code.git
  GIT_TAG        14.0.2  # Use the desired version tag
)
FetchContent_MakeAvailable(armadillo)

# Fetch LiquidDSP
FetchContent_Declare(
  liquid_dsp
  GIT_REPOSITORY https://github.com/jgaeddert/liquid-dsp.git
  GIT_TAG        v1.7.0
)
FetchContent_MakeAvailable(liquid_dsp)

# Set ROOT_DIR
get_filename_component(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}" ABSOLUTE)

# Include Armadillo
include_directories(SYSTEM "${armadillo_SOURCE_DIR}/include")

# Include LiquidDSP
include_directories(SYSTEM "${liquid_dsp_SOURCE_DIR}/include")

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
target_link_libraries(signalGen armadillo liquid)
