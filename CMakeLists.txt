cmake_minimum_required(VERSION 3.18)


project(IR_Library)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_subdirectory("IR_lib")
add_subdirectory("docs")
