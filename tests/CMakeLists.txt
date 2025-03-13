cmake_minimum_required(VERSION 3.10)

# Set compiler
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)

# Declare the project.
project(Tests)

# Set ROOT_DIR
get_filename_component(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

# Includes
include(FetchContent)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fetch Armadillo
FetchContent_Declare(
  armadillo
  GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code.git
  GIT_TAG        14.0.2  # Use the desired version tag
)
FetchContent_MakeAvailable(armadillo)

# Fetch Google Test
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.16.0
)
FetchContent_MakeAvailable(googletest)

# Fetch LiquidDSP
FetchContent_Declare(
  liquid_dsp
  GIT_REPOSITORY https://github.com/jgaeddert/liquid-dsp.git
  GIT_TAG        v1.7.0
)
FetchContent_MakeAvailable(liquid_dsp)

# Include Armadillo
include_directories(SYSTEM ${armadillo_SOURCE_DIR}/include)

# Include Google Test and Google Mock
include_directories(SYSTEM ${gtest_SOURCE_DIR}/include ${gmock_SOURCE_DIR}/include)

# Include LiquidDSP
include_directories(SYSTEM ${liquid_dsp_SOURCE_DIR}/include)

# Include project headers
include_directories(${ROOT_DIR}/inc)

# Collect all source files from src/ directory
file(GLOB CXX_SOURCES ${ROOT_DIR}/src/cpp/*.cpp)

# Collect all test files from tests/cpp/unit
file(GLOB UNIT_TEST_SOURCES ${ROOT_DIR}/tests/cpp/unit/*.cpp)

# Collect all test files from tests/cpp/integration
file(GLOB INTEGRATION_TEST_SOURCES ${ROOT_DIR}/tests/cpp/integration/*.cpp)

# Create the test executables
add_executable(unitTests ${UNIT_TEST_SOURCES} ${CXX_SOURCES})
add_executable(integrationTests ${INTEGRATION_TEST_SOURCES} ${CXX_SOURCES})

# Link with Armadillo, Google Test and pthread
target_link_libraries(unitTests armadillo gtest gtest_main gmock liquid pthread)
target_link_libraries(integrationTests armadillo gtest gtest_main gmock liquid pthread)

# Enable code coverage flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping -O0")

# Create report folder
file(MAKE_DIRECTORY ${ROOT_DIR}/htmlcov/cpp)

# Make coverage
add_custom_target(
    coverage
    COMMAND LLVM_PROFILE_FILE=coverageUnit.profraw ./unitTests
    COMMAND LLVM_PROFILE_FILE=coverageIntegration.profraw ./integrationTests
    COMMAND llvm-profdata merge -sparse coverageUnit.profraw coverageIntegration.profraw -o coverage.profdata
    COMMAND llvm-cov export ./unitTests ./integrationTests -instr-profile=coverage.profdata -format=lcov --ignore-filename-regex=".*armadillo_bits/.*" > ${ROOT_DIR}/htmlcov/cpp/coverage.lcov
    COMMAND genhtml ${ROOT_DIR}/htmlcov/cpp/coverage.lcov --output-directory ${ROOT_DIR}/htmlcov/cpp/
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
