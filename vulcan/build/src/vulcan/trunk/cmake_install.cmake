# Install script for directory: /media/skoll/FEM/vulcan/src/vulcan/trunk

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-m.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-m.O2")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-m.O2"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/bin/vulcan-bin/Vulcan-m.O2")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/bin/vulcan-bin" TYPE EXECUTABLE FILES "/media/skoll/FEM/vulcan/build/src/vulcan/trunk/Vulcan-m.O2")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-m.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-m.O2")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-m.O2")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-s.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-s.O2")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-s.O2"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/bin/vulcan-bin/Vulcan-s.O2")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/bin/vulcan-bin" TYPE EXECUTABLE FILES "/media/skoll/FEM/vulcan/build/src/vulcan/trunk/Vulcan-s.O2")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-s.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-s.O2")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-s.O2")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-t.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-t.O2")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-t.O2"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/bin/vulcan-bin/Vulcan-t.O2")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/bin/vulcan-bin" TYPE EXECUTABLE FILES "/media/skoll/FEM/vulcan/build/src/vulcan/trunk/Vulcan-t.O2")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-t.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-t.O2")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-t.O2")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tm.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tm.O2")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tm.O2"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/bin/vulcan-bin/Vulcan-tm.O2")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/bin/vulcan-bin" TYPE EXECUTABLE FILES "/media/skoll/FEM/vulcan/build/src/vulcan/trunk/Vulcan-tm.O2")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tm.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tm.O2")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tm.O2")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-ts.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-ts.O2")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-ts.O2"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/bin/vulcan-bin/Vulcan-ts.O2")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/bin/vulcan-bin" TYPE EXECUTABLE FILES "/media/skoll/FEM/vulcan/build/src/vulcan/trunk/Vulcan-ts.O2")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-ts.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-ts.O2")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-ts.O2")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tms.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tms.O2")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tms.O2"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/bin/vulcan-bin/Vulcan-tms.O2")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/bin/vulcan-bin" TYPE EXECUTABLE FILES "/media/skoll/FEM/vulcan/build/src/vulcan/trunk/Vulcan-tms.O2")
  if(EXISTS "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tms.O2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tms.O2")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/bin/vulcan-bin/Vulcan-tms.O2")
    endif()
  endif()
endif()

