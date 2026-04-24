# CMakeList.txt : CMake project for Thesis, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10)

project ("Thesis")

include(CTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

# Add source to this project's executable.
add_executable (Thesis "objects.cpp" "objects.h"
						"utils.cpp" "utils.h"
						"algorithms.cpp" "algorithms.h"
						"experiments.cpp" "experiments.h"
						"Thesis.cpp" "Thesis.h" "algorithmsNd.cpp" "algorithmsNd.h" "experimentsNd.cpp" "experimentsNd.h")

target_compile_options(Thesis PRIVATE /bigobj)

# add_executable(Example "example.cpp")

find_package(fmt CONFIG REQUIRED)
target_link_libraries(Thesis PRIVATE fmt::fmt)

find_package(Matplot++ CONFIG REQUIRED)
target_link_libraries(Thesis PRIVATE Matplot++::cimg  Matplot++::matplot)

include_directories("../lib/")
include_directories("../lib/eigen-5.0.0/")

# TODO: Add tests and install targets if needed.

# add_executable (Tests "objects.cpp" "objects.h" "objects_tests.cpp"
# 						"utils.cpp" "utils.h" "tests.cpp" "objects_tests.h")

add_test (NAME Tests COMMAND Tests)