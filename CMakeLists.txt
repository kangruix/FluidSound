cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(FluidSound LANGUAGES C CXX)

find_package(Eigen3 REQUIRED)

file(GLOB HEADER_FILES Source/*.h)
include_directories(Source/)

file(GLOB SRC_FILES Source/*.cpp Source/Integrators/*.cpp)

# Create executable only if FluidSound is being built as a standalone project 
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	add_executable(runFluidSound main.cpp ${SRC_FILES} ${HEADER_FILES})
	target_link_libraries(runFluidSound Eigen3::Eigen)

# Otherwise, build FluidSound as a static library
else ()
	add_library(${PROJECT_NAME} STATIC ${SRC_FILES} ${HEADER_FILES})
	target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

endif ()