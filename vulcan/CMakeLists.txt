#Make project file for Vulcan
cmake_minimum_required (VERSION 2.8)
project (Vulcan)
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

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # gfortran
    if(NOT restricted)
     set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fopenmp -cpp")
    else(NOT restricted)
     set (MAXELEM 3000)
     set (MAXNODE 1010)
     set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fopenmp -cpp -Drestricted -DMAXELEM=${MAXELEM} -DMAXNODE=${MAXNODE}")
    endif(NOT restricted)
    set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -cpp")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    # ifort
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -L$(MKLROOT)/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -openmp")
    set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set (EXECUTABLES "Vulcan-m.O2")
list (APPEND SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/scripts/vulcan")
list (APPEND SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/scripts/vulcan.py")

add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/vulcan/trunk")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/interpara/trunk")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/interfla/trunk")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/genren/trunk")
#add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/interexo/trunk")

if(NOT BUILD_MPACHECO)
  add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/intergplo/trunk")
  add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/intergplos/trunk")
else(NOT BUILD_MPACHECO)
  message("***********************************************")
  message("* SI QUIERE VOLVER A LA VERSION ORIGINAL DEBE *")
  message("* BORRAR LA CARPETA BUILD Y RECONSTRUIRLA     *")
  message("***********************************************")
  add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/intergplo/branches/mpacheco")
  add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/intergplos/branches/mpacheco")
endif(NOT BUILD_MPACHECO)
#add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/src/macout/trunk")

# install executables and scripts
set(CMAKE_INSTALL_PREFIX "/usr")
install (PROGRAMS ${SCRIPTS}
         DESTINATION "bin")

