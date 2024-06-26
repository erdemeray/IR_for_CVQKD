cmake_minimum_required(VERSION 3.18)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  set(CMAKE_CXX_FLAGS "/O2")
  set(CMAKE_CXX_FLAGS_RELEASE "/O2")
else()
  set(CMAKE_CXX_FLAGS "-O3")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

project(information_reconciliation)

include_directories("${CMAKE_SOURCE_DIR}/cpp/include/h_files")
include_directories("${CMAKE_SOURCE_DIR}/python")


set(PCM_DIR "PCM")
add_compile_definitions(PCM_DIR="${CMAKE_SOURCE_DIR}/cpp/${PCM_DIR}")


file (GLOB SOURCE_FILES "cpp/src/*.cpp")
file (GLOB HEADER_FILES "cpp/include/h_files/*.hpp")
file (GLOB PYTHON_FILES "python/*.cpp" "python/*.hpp")

# Set up such that XCode organizes the files
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )

find_package(Python COMPONENTS Interpreter Development REQUIRED)

#set(pybind11_DIR C:/Users/${USER}/anaconda3/Lib/site-packages/pybind11/share/cmake/pybind11) # Uncomment this line and change ${USER} if you have trouble with finding pybind11
find_package(pybind11 REQUIRED)

pybind11_add_module(information_reconciliation 
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)

target_link_libraries(information_reconciliation PUBLIC)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(information_reconciliation PRIVATE OpenMP::OpenMP_CXX)
endif()


install(TARGETS information_reconciliation
  COMPONENT python
  LIBRARY DESTINATION ${Python_SITELIB})

