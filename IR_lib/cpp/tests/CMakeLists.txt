cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(test_IR)

include_directories(/usr/local/include)
include_directories(../include)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../bin)

find_library(IR_LIB information_reconciliation HINTS ../build/) #

add_executable(test_IR src/use_case.cpp)

target_link_libraries(test_IR PUBLIC ${IR_LIB})


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_IR PRIVATE OpenMP::OpenMP_CXX)
endif()
