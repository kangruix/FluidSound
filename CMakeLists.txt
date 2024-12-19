cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)

project(FluidSound LANGUAGES C CXX)

find_package(Eigen3 3.4 REQUIRED)

file(GLOB H_FILES Source/*.h)
include_directories(Source/)

file(GLOB SRC_FILES Source/*.cpp Source/Integrators/*.cpp)

add_library(${PROJECT_NAME} STATIC ${SRC_FILES} ${H_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

# Create test executable only if FluidSound is being built as a standalone project 
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	add_executable(runFluidSound main.cpp)
	target_link_libraries(runFluidSound FluidSound Eigen3::Eigen)
endif ()