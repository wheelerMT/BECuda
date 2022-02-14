cmake_minimum_required(VERSION 3.21)
project(BECuda CUDA CXX)

set(CMAKE_CUDA_STANDARD 17)
find_package(CUDAToolkit REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)

set(INCLUDES include/BECuda.cuh include/constants.h include/grid.cuh include/wavefunction.cuh
             include/data.h include/kernels.cuh)
set(SOURCES src/grid.cu src/wavefunction.cu src/data.cpp src/kernels.cu)

add_subdirectory(HighFive)

add_library(BECuda STATIC ${INCLUDES} ${SOURCES} ${HIGHFIVEINCLUDE})
target_include_directories(${PROJECT_NAME} PUBLIC include/)
target_link_libraries(BECuda CUDA::cufft HighFive)

set_target_properties(BECuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

add_subdirectory(tests)
