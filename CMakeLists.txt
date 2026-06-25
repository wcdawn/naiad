cmake_minimum_required(VERSION 3.14)
project(naiad VERSION 0.1.0 DESCRIPTION "A one-dimensional SN neutron transport solver." LANGUAGES CXX)

add_executable(${PROJECT_NAME}.x)

target_sources(${PROJECT_NAME}.x
  PRIVATE
    src/analysis.cpp
    src/dense_matrix.cpp
    src/diffusion.cpp
    src/exception_handler.cpp
    src/geometry.cpp
    src/input.cpp
    src/linalg.cpp
    src/main.cpp
    src/math.cpp
    src/output.cpp
    src/quadrature_gauss_legendre.cpp
    src/quadrature_uniform.cpp
    src/quadrature1d.cpp
    src/result.cpp
    src/timer.cpp
    src/transport.cpp
    src/writer.cpp
    src/xslibrary.cpp
)

set_property(TARGET ${PROJECT_NAME}.x PROPERTY CXX_STANDARD 23)

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME}.x OpenMP::OpenMP_CXX)

target_compile_options(${PROJECT_NAME}.x PRIVATE -Wall -Wextra -O3)
