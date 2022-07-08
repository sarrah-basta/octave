cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(example)

include_directories("/usr/local/include" "/usr/include/octave-8.0.0/octave")
link_directories("/usr/local/lib" "/usr/lib/octave/8.0.0")
add_executable(app example.cpp nvector_col.cpp nvector_col.h)
target_link_libraries(app sundials_generic octave)
