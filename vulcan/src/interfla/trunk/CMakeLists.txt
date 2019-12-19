cmake_minimum_required (VERSION 2.8)
project (Interfla)
enable_language (Fortran)

add_executable ("interfla" "interfla.f")

install(TARGETS "interfla" DESTINATION "bin")
