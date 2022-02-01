cmake_minimum_required(VERSION 3.21)
project(CodeTesting CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(${PROJECT_NAME} tests.cu)

target_link_libraries(${PROJECT_NAME} PUBLIC BECuda)
