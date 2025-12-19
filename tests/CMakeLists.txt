cmake_minimum_required(VERSION 3.10)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Declare the project.
project(Tests)

# Set ROOT_DIR
get_filename_component(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

# Includes
include(FetchContent)

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

# Targets
set(INCLUDE_DIR "${ROOT_DIR}/inc")
file(GLOB CXX_SOURCES "${ROOT_DIR}/src/cpp/*.cpp")
list(REMOVE_ITEM CXX_SOURCES "${ROOT_DIR}/src/cpp/main.cpp")

set(TESTS_DIR "${ROOT_DIR}/tests/cpp")

# Unit
file(GLOB UNIT_TEST_SOURCES "${TESTS_DIR}/unit/*.cpp")
add_executable(unitTests "${CXX_SOURCES}" "${UNIT_TEST_SOURCES}")

target_include_directories(
  unitTests PRIVATE
  "${INCLUDE_DIR}"
  "${armadillo_SOURCE_DIR}/include"
  "${gtest_SOURCE_DIR}/include"
  "${gmock_SOURCE_DIR}/include"
  "${liquid_dsp_SOURCE_DIR}/include"
)

target_link_libraries(unitTests PRIVATE armadillo gtest gtest_main gmock liquid pthread)

add_custom_target(
  unit
  COMMAND "./unitTests"
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# Integration
file(GLOB INTEGRATION_TEST_SOURCES "${ROOT_DIR}/tests/cpp/integration/*.cpp")
add_executable(integrationTests "${INTEGRATION_TEST_SOURCES}" "${CXX_SOURCES}")

target_include_directories(
  integrationTests PRIVATE
  "${INCLUDE_DIR}"
  "${armadillo_SOURCE_DIR}/include"
  "${gtest_SOURCE_DIR}/include"
  "${gmock_SOURCE_DIR}/include"
  "${liquid_dsp_SOURCE_DIR}/include"
)

target_link_libraries(integrationTests PRIVATE armadillo gtest gtest_main gmock liquid pthread)

add_custom_target(
  integration
  COMMAND "./integrationTests"
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# Enable code coverage flags
target_compile_options(unitTests PRIVATE --coverage -O0 -g)
target_link_options(unitTests PRIVATE --coverage)

target_compile_options(integrationTests PRIVATE --coverage -O0 -g)
target_link_options(integrationTests PRIVATE --coverage)

# Create report folder
set(REPORT_DIR "${ROOT_DIR}/reports/htmlcov/cpp")
file(MAKE_DIRECTORY "${REPORT_DIR}")

# Make coverage
add_custom_target(
  coverage
  COMMAND ./unitTests
  COMMAND ./integrationTests
  COMMAND lcov --capture --directory . -o raw.info --ignore-errors inconsistent,inconsistent --ignore-errors mismatch,mismatch --rc geninfo_unexecuted_blocks=1 > /dev/null
  COMMAND lcov --extract raw.info "${ROOT_DIR}/inc/*" "${ROOT_DIR}/src/cpp/*" -o coverage.info --ignore-errors unused,unused > /dev/null
  COMMAND sed -i "/D0Ev/d" coverage.info
  COMMAND genhtml coverage.info -o "${REPORT_DIR}"
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
)
