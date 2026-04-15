cmake_minimum_required(VERSION 3.28)

# ===== Project =====
project(ProgrammableSignalGenerator)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ===== Root Directory =====
get_filename_component(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}" ABSOLUTE)

# ===== Set Variables =====
set(APP_DIR "${ROOT_DIR}/app")
set(OUTPUT_DIR "${ROOT_DIR}/bin" CACHE PATH "Output Directory")
set(MODULE_NAME "signalgen" CACHE STRING "Output Module Name")
set(INCLUDE_DIR "${ROOT_DIR}/inc" CACHE PATH "C++ Headers")
set(SOURCES_DIR "${ROOT_DIR}/src/cpp" CACHE PATH "C++ Sources")
file(GLOB CXX_SOURCES "${SOURCES_DIR}/*.cpp")

# ===== Additional Dependencies =====
find_package(armadillo REQUIRED)
find_package(liquid-dsp REQUIRED)

# ===== Add Additional Sources =====
list(APPEND CXX_SOURCES "${APP_DIR}/main.cpp")

# ===== Executable =====
add_executable(signalgen_exe "${CXX_SOURCES}")

target_include_directories(signalgen_exe PRIVATE "${INCLUDE_DIR}")

target_link_libraries(signalgen_exe armadillo::armadillo liquid-dsp::liquid-dsp)

set_target_properties(
  signalgen_exe PROPERTIES
  OUTPUT_NAME "${MODULE_NAME}"
  RUNTIME_OUTPUT_DIRECTORY "${OUTPUT_DIR}"
)
