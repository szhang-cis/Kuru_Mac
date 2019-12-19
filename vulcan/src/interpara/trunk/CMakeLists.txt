cmake_minimum_required (VERSION 2.8)
project (Interpara)
enable_language (Fortran)

add_executable ("interpara_for" "interpara.for")

install(TARGETS "interpara_for" DESTINATION "/usr/bin/vulcan-bin")
install(PROGRAMS "interpara" DESTINATION "bin")
