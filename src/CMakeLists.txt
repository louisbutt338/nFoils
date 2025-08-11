cmake_minimum_required(VERSION 3.5 FATAL_ERROR)               
project(nfoils CXX)
add_library(functions MODULE functions.cpp)
target_link_libraries(functions)