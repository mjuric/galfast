cmake_minimum_required(VERSION 2.6)
project(galfast)

#
# git-based automatic versioning
#
add_custom_target(
  version-gen
  COMMAND "${PROJECT_SOURCE_DIR}/git-version-gen" "${PROJECT_BINARY_DIR}/src/galfast_version.h" "${PROJECT_SOURCE_DIR}/.version"
  COMMENT "Generating version information"
)
execute_process(
  COMMAND "${PROJECT_SOURCE_DIR}/git-version-gen" "${PROJECT_BINARY_DIR}/src/galfast_version.h" "${PROJECT_SOURCE_DIR}/.version"
)

#
# Check requirements
#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

### fmemopen
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists(fmemopen HAVE_FMEMOPEN)

### libpeyton
find_package(libpeyton REQUIRED)
include_directories(${libpeyton_INCLUDE_DIR})

### GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

### Boost
find_package(Boost REQUIRED COMPONENTS regex iostreams thread)
include_directories(${Boost_INCLUDE_DIRS})

### CFITSIO
find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIR})

### ZLIB
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR})

### ZLIB
find_package(BZip2 REQUIRED)
include_directories(${BZip2_INCLUDE_DIR})

### CUDA
find_package(CUDA 2.3 REQUIRED)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_13)

if( NOT "x${GCC_ROOT}" STREQUAL "x" )
	set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-ccbin=${GCC_ROOT}/bin")
	message(STATUS "Note: CUDA will use host compiler from ${GCC_ROOT}")
endif( NOT "x${GCC_ROOT}" STREQUAL "x" )

if( "x${CUDA_BUILD_EMULATION}" STREQUAL xON )
	set(CUDA_DEVEMU 1)
	message(STATUS "Note: Targeting CUDA code for device emulation")
endif( "x${CUDA_BUILD_EMULATION}" STREQUAL xON )

#
# Data Directories
#
set(DATADIRSUFFIX "share/${CMAKE_PROJECT_NAME}")
set(DATADIR "${CMAKE_INSTALL_PREFIX}/${DATADIRSUFFIX}")

install(FILES
	data/readme.txt

	data/safeprimes32.txt

	data/fbMr.txt

	data/SDSSugriz.photosys.txt
	data/SDSSugriz.SDSSg.photoerr.txt data/SDSSugriz.SDSSi.photoerr.txt data/SDSSugriz.SDSSr.photoerr.txt data/SDSSugriz.SDSSu.photoerr.txt data/SDSSugriz.SDSSz.photoerr.txt

	data/LSSTugrizy.photosys.txt
	data/LSSTugrizy.LSSTi.photoerr.txt data/LSSTugrizy.LSSTr.photoerr.txt data/LSSTugrizy.LSSTu.photoerr.txt data/LSSTugrizy.LSSTz.photoerr.txt data/LSSTugrizy.LSSTy.photoerr.txt data/LSSTugrizy.LSSTg.photoerr.txt
DESTINATION ${DATADIRSUFFIX})

#
# configure a header file to pass some of the CMake settings to the source code
#
configure_file(
  "${PROJECT_SOURCE_DIR}/galfast_config.h.in"
  "${PROJECT_BINARY_DIR}/src/galfast_config.h"
  ESCAPE_QUOTES
)
include_directories("${PROJECT_BINARY_DIR}/src")

#
# galfast sources
#
include_directories("${PROJECT_SOURCE_DIR}/src/common")
include_directories("${PROJECT_SOURCE_DIR}/src/modules")
include_directories("${PROJECT_SOURCE_DIR}/src/skygen")
cuda_add_executable(galfast

  src/galfast.cpp
  src/pipeline.cpp
  src/kernels_cpu.cpp
  src/footprint.cpp
  src/tests.cpp

  src/kernels_gpu.cu

  src/gpc/gpc.c

  src/skygen/model_expDisk.cpp
  src/skygen/model_powerLawEllipsoid.cpp
  src/skygen/os_skygen.cpp
  src/skygen/model_brokenPowerLaw.cpp
  src/skygen/model_densityCube.cpp
  src/skygen/model_J08.cpp

  src/common/gpc_cpp.cpp
  src/common/spline.cpp
  src/common/fitter.cpp
  src/common/transform.cpp
  src/common/otable.cpp
  src/common/cux.cpp
  src/common/binarystream.cpp
  src/common/io.cpp

  src/modules/Bond2010_host.cpp
  src/modules/photometricErrors_host.cpp
  src/modules/GaussianFeH_host.cpp
  src/modules/gal2other_host.cpp
  src/modules/fixedFeH_host.cpp
  src/modules/unresolvedMultiples_host.cpp
  src/modules/kinTMIII_host.cpp
  src/modules/photometry_host.cpp
  src/modules/vel2pm_host.cpp
  src/modules/FeH_host.cpp

)

target_link_libraries(galfast ${Boost_LIBRARIES})
target_link_libraries(galfast ${CFITSIO_LIBRARIES})
target_link_libraries(galfast ${GSL_LIBRARIES})
target_link_libraries(galfast ${libpeyton_LIBRARIES})
target_link_libraries(galfast ${ZLIB_LIBRARIES})
target_link_libraries(galfast ${BZip2_LIBRARIES})

add_dependencies(galfast version-gen)

#
# install info
#
install (TARGETS galfast DESTINATION bin)

#
# package creation info
#   - this generates a "dist" target to create binary and source packages
#

execute_process(
  COMMAND "${PROJECT_SOURCE_DIR}/git-untracked-files" "${PROJECT_SOURCE_DIR}"
  OUTPUT_VARIABLE CPACK_IGNORE_FILES
)

add_custom_target(dist
	rm -f CPackConfig.make CPackSourceConfig.make &&
	cmake "${CMAKE_CURRENT_SOURCE_DIR}" &&
	make package_source &&
	make package
)

file(READ "${PROJECT_SOURCE_DIR}/.version" VLIST)
STRING(REGEX MATCHALL "[^ ]+" VLIST ${VLIST})
list(GET VLIST 0 CPACK_PACKAGE_FILE_NAME)
list(GET VLIST 1 CPACK_PACKAGE_VERSION_MAJOR)
list(GET VLIST 2 CPACK_PACKAGE_VERSION_MINOR)
list(GET VLIST 3 CPACK_PACKAGE_VERSION_PATCH)
list(GET VLIST 4 CPACK_RPM_PACKAGE_RELEASE)
file(READ "${PROJECT_SOURCE_DIR}/.version.fn" FULL_VERSION)

include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set (CPACK_GENERATOR "TBZ2")
set (CPACK_SOURCE_GENERATOR "TBZ2")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "galfast - fast mock catalog generator")
set (CPACK_PACKAGE_DEFAULT_LOCATION "/usr")
set (CPACK_PACKAGE_VENDOR "Mario Juric")
set (CPACK_PACKAGE_CONTACT "Mario Juric <mjuric@ias.edu>")
set (CPACK_PACKAGE_FILE_NAME "galfast-${FULL_VERSION}.bin")
set (CPACK_SOURCE_PACKAGE_FILE_NAME "galfast-${FULL_VERSION}")
set (CPACK_SOURCE_IGNORE_FILES ${CPACK_IGNORE_FILES})
include (CPack)
