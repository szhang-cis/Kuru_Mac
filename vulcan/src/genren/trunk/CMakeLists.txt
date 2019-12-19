cmake_minimum_required (VERSION 2.8)
project (Genren)
enable_language (Fortran)

add_executable ("genren.out" "genren.f")
add_executable ("malla3d" "malla3d.f")
add_executable ("mlit-mesh" "mlit-mesh.f")
add_executable ("9nto4n" "9nto4n.f")

install(TARGETS "genren.out" DESTINATION "/usr/bin/vulcan-bin")
install(PROGRAMS "genren" DESTINATION "bin")
install(TARGETS "malla3d" DESTINATION "bin")
install(TARGETS "mlit-mesh" DESTINATION "bin")
install(TARGETS "9nto4n" DESTINATION "bin")
