cmake_minimum_required(VERSION 3.18)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(information_reconciliation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set the output directory for the build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories("${CMAKE_SOURCE_DIR}/cpp/include/h_files")
include_directories("${CMAKE_SOURCE_DIR}/python")


set(PCM_DIR "PCM")
add_compile_definitions(PCM_DIR="${CMAKE_SOURCE_DIR}/cpp/${PCM_DIR}")

set(PCM_NAME "R_0p2_R_0p01_RA_SPPCOM_2")
add_compile_definitions(PCM_NAME="R_0p2_R_0p01_RA_SPPCOM_2")

# Detect conda environment - COMMENT IF YOU DON'T USE CONDA
if(DEFINED ENV{CONDA_PREFIX})
    set(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})
    message(STATUS "Conda environment detected: $ENV{CONDA_PREFIX}")
else()
    message(WARNING "Conda environment not detected. You may need to activate your environment.")
endif()

# Set Python variables
set(Python_ROOT_DIR ${CMAKE_PREFIX_PATH})
set(Python_USE_STATIC_LIBS FALSE)
set(Python_FIND_STRATEGY LOCATION)
set(Python_FIND_REGISTRY NEVER)
set(Python_FIND_FRAMEWORK NEVER)

# Find Python package
find_package(Python 3.9 EXACT REQUIRED COMPONENTS Interpreter Development)

# Ensure we're using the correct Python version
if(${Python_VERSION} VERSION_LESS 3.9)
  message(FATAL_ERROR "Python version 3.9 or newer is required. Found: ${Python_VERSION}")
endif()

# Print Python variables for debugging
message(STATUS "Python_FOUND: ${Python_FOUND}")
message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")
message(STATUS "Python_INCLUDE_DIRS: ${Python_INCLUDE_DIRS}")
message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES}")

# Ensure pybind11 uses the correct Python version
set(PYBIND11_PYTHON_VERSION 3.9)
set(PYBIND11_PYTHON_EXECUTABLE ${Python_EXECUTABLE})
set(PYBIND11_PYTHON_INCLUDE_DIR ${Python_INCLUDE_DIR})
set(PYBIND11_PYTHON_LIBRARIES ${Python_LIBRARIES})
set(PYBIND11_FIND_PYTHON ON)

# Include pybind11
include(FetchContent)

FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/v2.13.0.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

set(PYBIND11_NOPYTHON FALSE)
set(PYBIND11_PYTHON_VERSION 3.9)

# Pass the already found Python variables to pybind11
set(PYTHON_INCLUDE_DIR ${Python_INCLUDE_DIRS})
set(PYTHON_LIBRARY ${Python_LIBRARIES})
set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})

FetchContent_MakeAvailable(pybind11)

# Gather source files
file(GLOB SOURCE_FILES "cpp/src/*.cpp")
file(GLOB HEADER_FILES "cpp/include/h_files/*.hpp")
file(GLOB PYTHON_FILES "python/*.cpp" "python/*.hpp")

# Create the pybind11 module
pybind11_add_module(information_reconciliation 
    ${SOURCE_FILES}
    ${HEADER_FILES}
    ${PYTHON_FILES}
)

# Include directories
target_include_directories(information_reconciliation PRIVATE 
    ${CMAKE_SOURCE_DIR}/cpp/include/h_files
    ${CMAKE_SOURCE_DIR}/python
    ${Python_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(information_reconciliation PRIVATE ${Python_LIBRARIES})


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(information_reconciliation PRIVATE OpenMP::OpenMP_CXX)
endif()

# Install the module
install(TARGETS information_reconciliation
    COMPONENT python
    LIBRARY DESTINATION "${Python_SITEARCH}"
    RUNTIME DESTINATION "${Python_SITEARCH}"
)
