cmake_minimum_required(VERSION 3.28)

# ===== Project =====
project(Tests)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ===== Root Directory =====
get_filename_component(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

# ===== Set Variables =====
set(INCLUDE_DIR "${ROOT_DIR}/inc" CACHE PATH "C++ Headers")
set(SOURCES_DIR "${ROOT_DIR}/src/cpp" CACHE PATH "C++ Sources")
file(GLOB CXX_SOURCES "${SOURCES_DIR}/*.cpp")

# ===== Additional Dependencies =====
find_package(armadillo REQUIRED)
find_package(GTest REQUIRED)
find_package(liquid-dsp REQUIRED)

# ===== Targets =====
set(TESTS_DIR "${ROOT_DIR}/tests/cpp")

# ----- Unit Tests -----
file(GLOB UNIT_TEST_SOURCES "${TESTS_DIR}/unit/*.cpp")
add_executable(unitTests "${CXX_SOURCES}" "${UNIT_TEST_SOURCES}")

target_include_directories(unitTests PRIVATE "${INCLUDE_DIR}")

target_link_libraries(unitTests PRIVATE
  armadillo::armadillo
  GTest::gtest
  GTest::gtest_main
  GTest::gmock
  liquid-dsp::liquid-dsp
  pthread
)

add_custom_target(
  unit
  COMMAND "./unitTests"
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# ----- Integration Tests -----
file(GLOB INTEGRATION_TEST_SOURCES "${ROOT_DIR}/tests/cpp/integration/*.cpp")
add_executable(integrationTests "${INTEGRATION_TEST_SOURCES}" "${CXX_SOURCES}")

target_include_directories(integrationTests PRIVATE "${INCLUDE_DIR}")

target_link_libraries(integrationTests PRIVATE
  armadillo::armadillo
  GTest::gtest
  GTest::gtest_main
  GTest::gmock
  liquid-dsp::liquid-dsp
  pthread
)

add_custom_target(
  integration
  COMMAND "./integrationTests"
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# ===== Coverage =====
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
