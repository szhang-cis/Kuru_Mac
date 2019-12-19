cmake_minimum_required (VERSION 2.8)
project (Interexo)
enable_language (Fortran)

add_executable ("interexo" interexo.f)

target_link_libraries("interexo" vulexo2for.so vulexo.so)

install(TARGETS "interexo" DESTINATION "/usr/bin/")
