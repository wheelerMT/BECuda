cmake_minimum_required(VERSION 3.21)
project(BECuda CUDA)

set(CMAKE_CUDA_STANDARD 17)
find_package(CUDAToolkit REQUIRED)

set(INCLUDES include/BECuda.cuh include/constants.h include/grid.cuh include/wavefunction.cuh)
set(SOURCES src/grid.cu src/wavefunction.cu)

add_library(BECuda STATIC ${INCLUDES} ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC include/)
target_link_libraries(BECuda CUDA::cufft)

set_target_properties(BECuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

add_subdirectory(tests)
