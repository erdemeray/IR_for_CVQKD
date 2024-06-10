cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(information_reconciliation VERSION 0.1.0)

# Include dir
include_directories(/usr/local/include)

# Src
AUX_SOURCE_DIRECTORY(src SRC_FILES)

# Directory Definitions
set(PROJECT_SOURCE_DIR "src")
set(PROJECT_INCLUDE_DIR "include/h_files")
set(PCM_DIR "PCM")

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Source files
set(SOURCE_FILES
    ${PROJECT_INCLUDE_DIR}/reconciliation.hpp
    ${PROJECT_INCLUDE_DIR}/utilities.hpp	
    ${PROJECT_SOURCE_DIR}/reconciliation.cpp
    ${PROJECT_SOURCE_DIR}/utilities.cpp
)

# Add compile definitions for PCM locations
add_compile_definitions(PCM_DIR="${CMAKE_SOURCE_DIR}/${PCM_DIR}")

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES})

# Add library
add_library(information_reconciliation SHARED ${SOURCE_FILES})

# Include directories
target_include_directories(information_reconciliation PRIVATE include/information_reconciliation)

# Install
install(TARGETS information_reconciliation DESTINATION lib)

# Install the headers
install(FILES include/information_reconciliation DESTINATION include)

# Create base directory
install(DIRECTORY include/h_files DESTINATION include)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(information_reconciliation PRIVATE OpenMP::OpenMP_CXX)
endif()





