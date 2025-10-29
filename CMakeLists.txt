# CMakeList.txt : CMake project for Thesis, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10)

project ("Thesis")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add source to this project's executable.
add_executable (Thesis "Thesis.cpp" "Thesis.h")
add_executable(Example "example.cpp")

find_package(fmt CONFIG REQUIRED)
target_link_libraries(Thesis PRIVATE fmt::fmt)

find_package(Matplot++ CONFIG REQUIRED)
target_link_libraries(Thesis PRIVATE Matplot++::cimg  Matplot++::matplot)

include_directories("../lib/")

# TODO: Add tests and install targets if needed.