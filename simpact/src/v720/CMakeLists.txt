#Make project file for Vulcan
cmake_minimum_required (VERSION 2.8)
project (simpact)
enable_language (Fortran)

find_package(OpenMP)
if(OPENMP_FOUND)
  set (CMAKE_F_FLAGS "${CMAKE_F_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
         "Choose the type of build, options are: None Debug Release."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set (CMAKE_Fortran_FLAGS_RELEASE "-O1 -fopenmp -cpp")
#set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fopenmp -cpp")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inclu)
set(CMAKE_INSTALL_PREFIX "/usr/local/")


file(GLOB_RECURSE simpact_SOURCES "auxil/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "beame/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "cont4/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "drawb/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "input/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "lbst/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "loads/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "mainp/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "matdb/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "nbst/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "remesh/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "rigid/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "scpr03/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "solag/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "solve/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "spot/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "ssolag/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})
file(GLOB_RECURSE simpact_SOURCES "truss/*.f90")
list (APPEND OBJECT_LIST ${GLOB_RECURSE})

#list(REMOVE_DUPLICATES ${OBJECT_LIST})

add_executable ("simpact.o" ${OBJECT_LIST})
