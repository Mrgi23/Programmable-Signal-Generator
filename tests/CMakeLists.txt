cmake_minimum_required(VERSION 3.10)

# Name of the project
project(Tests)

# Set ROOT_DIR
get_filename_component(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

# Includes
include(FetchContent)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fetch Google Test if not installed
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)

# Prevent overriding parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Fetch and add Google Test
FetchContent_MakeAvailable(googletest)

# Include Google Test and Google Mock
include_directories(${gtest_SOURCE_DIR}/include ${gmock_SOURCE_DIR}/include)

# Set compiler
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)

# Include project headers
include_directories(${ROOT_DIR}/inc)

# Collect all source files from src/ directory
file(GLOB CXX_SOURCES ${ROOT_DIR}/src/cpp/*.cpp)

# Collect all test files from tests/c++/
file(GLOB TEST_SOURCES ${ROOT_DIR}/tests/cpp/*.cpp)

# Create a single test executable
add_executable(runAllTests.exe ${TEST_SOURCES} ${CXX_SOURCES})

# Link with Google Test and pthread
target_link_libraries(runAllTests.exe gtest gtest_main gmock pthread)

# Enable code coverage flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping -O0")

# Create report folder
file(MAKE_DIRECTORY ${ROOT_DIR}/htmlcov/cpp)

# Make coverage
add_custom_target(coverage
    COMMAND LLVM_PROFILE_FILE=coverage.profraw ./runAllTests.exe
    COMMAND llvm-profdata merge -sparse coverage.profraw -o coverage.profdata
    COMMAND llvm-cov export ./runAllTests.exe -instr-profile=coverage.profdata -format=lcov > ${ROOT_DIR}/htmlcov/cpp/coverage.lcov
    COMMAND genhtml ${ROOT_DIR}/htmlcov/cpp/coverage.lcov --output-directory ${ROOT_DIR}/htmlcov/cpp/
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
