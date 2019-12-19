# CMake project file for Simpact

##################################################
# Define the project and the depencies that it has
##################################################
cmake_minimum_required (VERSION 2.8)
project (simpact)
enable_language (Fortran)

# Set de Simpact version
#set (VERSION 1.0.0)

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
set (curvasexe curvas)
#set (stpexe stp)

# Define the library name
set (stampactlib stampact)

# Define some directories
set (src ${CMAKE_SOURCE_DIR}/src)
set (lib ${CMAKE_SOURCE_DIR}/lib)
set (bin ${CMAKE_SOURCE_DIR}/bin)
set (src_curvas ${src}/curvas)
#set (src_stp ${src}/stp)
#set (src_v720 ${src}/v720)
set (src_stampact ${src}/stampact)

# Have the .mod files placed in the lib folder
set (CMAKE_Fortran_MODULE_DIRECTORY ${lib})

# The source for the stampact library and have it placed in the lib folder
add_subdirectory (${src_stampact} ${lib})

# The source for the simpact binary and have it placed in the bin folder
add_subdirectory (${src_curvas} ${bin})
#add_subdirectory (${src_stp} ${bin})
