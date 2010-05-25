# - Try to find libpeyton
# Once done this will define
#
#  libpeyton_FOUND - system has LIBPEYTON
#  libpeyton_INCLUDE_DIR - the LIBPEYTON include directory
#  libpeyton_LIBRARIES - Link these to use LIBPEYTON
#
#  If libpeyton_ROOT is set, look for it there
#

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

  find_path(libpeyton_INCLUDE_DIR astro/peyton_version.h
    PATHS ${libpeyton_ROOT}/include
  )

  find_library(libpeyton_LIBRARIES NAMES peyton
    PATHS ${libpeyton_ROOT}/lib
  )

  if(libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)
    set(libpeyton_FOUND TRUE)
  else (libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)
    set(libpeyton_FOUND FALSE)
  endif(libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)

  if (libpeyton_FOUND)
    if (NOT libpeyton_FIND_QUIETLY)
      message(STATUS "Found libpeyton ${libpeyton_VERSION_STRING}: ${libpeyton_LIBRARIES}")
    endif (NOT libpeyton_FIND_QUIETLY)
  else (libpeyton_FOUND)
    if (libpeyton_FIND_REQUIRED)
      message(SEND_ERROR "libpeyton not found.")
    else (libpeyton_FIND_REQUIRED)
      message(STATUS "libpeyton not found.")
    endif (libpeyton_FIND_REQUIRED)
  endif (libpeyton_FOUND)

  mark_as_advanced(libpeyton_INCLUDE_DIR libpeyton_LIBRARIES)

endif (libpeyton_INCLUDE_DIR AND libpeyton_LIBRARIES)
