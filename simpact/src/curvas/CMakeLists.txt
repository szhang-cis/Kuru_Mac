########################################
# Set up how to compile the source files
########################################

# Add the source files
set (curvas_src ${src_curvas}/history.f90
                ${src_curvas}/histd.f90
                ${src_curvas}/stacpvt.f90
                ${src_curvas}/stadisp.f90
                ${src_curvas}/stastre.f90
                ${src_curvas}/statcfr.f90
)

# Define the executable in terms of the source files
add_executable (${curvasexe} ${curvas_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links simpact to the simpact library
target_link_libraries (${curvasexe} ${stampactlib})
