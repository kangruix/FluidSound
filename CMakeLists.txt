cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(FluidSound LANGUAGES C CXX)

option(FORCE_CPU "Force CPU-only compilation (even if GPU is available)" ON)

find_package(Eigen3 REQUIRED)
find_package(CUDA)

file(GLOB HEADER_FILES Source/*.h)
include_directories(Source/)

if (CUDA_FOUND AND NOT FORCE_CPU)
	message(STATUS "Building for GPU")
	enable_language(CUDA)

	file(GLOB SRC_FILES Source/*.cpp Source/*.cu)
	add_executable(${PROJECT_NAME} main.cpp ${SRC_FILES} ${HEADER_FILES})

	target_compile_definitions(${PROJECT_NAME} PRIVATE USE_CUDA)
	target_link_libraries(${PROJECT_NAME} Eigen3::Eigen cublas)

else ()
	message(STATUS "Building for CPU")

	file(GLOB SRC_FILES Source/*.cpp)
	add_executable(${PROJECT_NAME} main.cpp ${SRC_FILES} ${HEADER_FILES})

	target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
endif ()