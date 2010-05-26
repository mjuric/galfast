# Try to find libpeyton. Look for it in
#  1) Standard include/library paths
#  2) libpeyton_ROOT/lib and libpeyton_ROOT/include
#  3) libpeyton_ROOT (for lib) and libpeyton_ROOT/include plus libpeyton_ROOT/../include
#     (this is the usual setup if libpeyton_ROOT is an in-source build directory)
#
# Once done this will define
#
#  libpeyton_FOUND - system has LIBPEYTON
#  libpeyton_INCLUDE_DIR - the LIBPEYTON include directory
#  libpeyton_LIBRARIES - Link these to use LIBPEYTON
#
#  If libpeyton_ROOT is set, look for it there
#

# Copyright (c) 2010, Mario Juric <mjuric@ias.edu>
# Copyright (c) 2006, Jasem Mutlaq <mutlaqja@ikarustech.com>
# Based on FindLibfacile by Carsten Niehaus, <cniehaus@gmx.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if (libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)

  # in cache already
  set(libpeyton_FOUND TRUE)
  message(STATUS "Found libpeyton: ${libpeyton_LIBRARIES}")

else (libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)

    find_path(libpeyton_INCLUDE_DIR astro/version.h
      PATHS ${libpeyton_ROOT}/include
    )

    find_library(libpeyton_LIBRARIES NAMES peyton
      PATHS ${libpeyton_ROOT}/lib
    )

  if(NOT (libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES) AND libpeyton_ROOT)
    # Reinterpret libpeyton_ROOT as referring to a build directory
    # within the root source directory, and try to find include+lib files

    # Find generated include files
    find_path(libpeyton_INCLUDE_DIR_static astro/peyton_version.h
      PATHS ${libpeyton_ROOT}/include
    )

    # Find static include files
    find_path(libpeyton_INCLUDE_DIR_gen astro/version.h
      PATHS ${libpeyton_ROOT}/../include
    )

    # Set libpeyton_INCLUDE_DIR only if both static and generated include
    # files were found.
    if (libpeyton_INCLUDE_DIR_gen AND libpeyton_INCLUDE_DIR_static)
      set (libpeyton_INCLUDE_DIR ${libpeyton_INCLUDE_DIR_static} ${libpeyton_INCLUDE_DIR_gen})
    endif (libpeyton_INCLUDE_DIR_gen AND libpeyton_INCLUDE_DIR_static)

    # Find the library
    find_library(libpeyton_LIBRARIES NAMES peyton
      PATHS ${libpeyton_ROOT}
    )
  endif(NOT (libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES) AND libpeyton_ROOT)


  if(libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)
    set(libpeyton_FOUND TRUE)
  else (libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)
    set(libpeyton_FOUND FALSE)
  endif(libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)

  if (libpeyton_FOUND)
    if (NOT libpeyton_FIND_QUIETLY)
      message(STATUS "Found libpeyton ${libpeyton_VERSION_STRING}: ${libpeyton_INCLUDE_DIR}")
      message(STATUS "Found libpeyton ${libpeyton_VERSION_STRING}: ${libpeyton_LIBRARIES}")
    endif (NOT libpeyton_FIND_QUIETLY)
  else (libpeyton_FOUND)
    if (libpeyton_FIND_REQUIRED)
      message(FATAL_ERROR "libpeyton not found.")
    else (libpeyton_FIND_REQUIRED)
      message(STATUS "libpeyton not found.")
    endif (libpeyton_FIND_REQUIRED)
  endif (libpeyton_FOUND)

  mark_as_advanced(libpeyton_INCLUDE_DIR libpeyton_LIBRARIES)

endif (libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)
