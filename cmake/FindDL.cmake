# - Try to find libdl
# Once done this will define
#  DL_FOUND - System has libdl
#  DL_INCLUDE_DIRS - The libdl include directories
#  DL_LIBRARY_DIRS - The library directories needed to use libdl
#  DL_LIBRARIES    - The libraries needed to use libdl

if (DL_INCLUDES AND DL_LIBRARIES)
  set(DL_FIND_QUIETLY TRUE)
endif (DL_INCLUDES AND DL_LIBRARIES)


FIND_PACKAGE(PkgConfig)
IF( PKG_CONFIG_FOUND )

  pkg_check_modules( DL QUIET dl)

ENDIF( PKG_CONFIG_FOUND )


IF( DL_FOUND )
  IF(NOT DL_FIND_QUIETLY)
    MESSAGE(STATUS "package dl found")
  ENDIF()
ELSE( DL_FOUND )
  FIND_PATH(DL_INCLUDE_DIRS NAMES dlfcn.h
            HINTS  $ENV{DL_HOME}
                   "/usr/include/" )
           
  FIND_LIBRARY(DL_LIBRARY 
               dl ltdl
                HINTS "$ENV{DL_HOME}/lib"
                      "/usr/lib" )

  SET(DL_INCLUDE_DIRS ${DL_INCLUDE_DIR})
  SET(DL_LIBRARIES ${DL_LIBRARY})

  INCLUDE(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set LIBDL_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(DL  DEFAULT_MSG
                                    DL_INCLUDE_DIRS
                                    DL_LIBRARIES)

  MARK_AS_ADVANCED(DL_INCLUDE_DIRS DL_LIBRARIES)
  
  IF( DL_FOUND AND NOT DL_FIND_QUIETLY )
    MESSAGE(STATUS "package dl found")
  ENDIF()
ENDIF( DL_FOUND )








# - Find libdl
# Find the native LIBDL includes and library
#
# LIBDL_INCLUDE_DIR - where to find dlfcn.h, etc.
# LIBDL_LIBRARIES - List of libraries when using libdl.
# LIBDL_FOUND - True if libdl found.


IF (LIBDL_INCLUDE_DIR)
# Already in cache, be silent
SET(LIBDL_FIND_QUIETLY TRUE)
ENDIF (LIBDL_INCLUDE_DIR)

FIND_PATH(LIBDL_INCLUDE_DIR dlfcn.h)

SET(LIBDL_NAMES dl libdl ltdl libltdl)
FIND_LIBRARY(LIBDL_LIBRARY NAMES ${LIBDL_NAMES} )

# handle the QUIETLY and REQUIRED arguments and set LIBDL_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LibDL DEFAULT_MSG LIBDL_LIBRARY LIBDL_INCLUDE_DIR)

IF(LIBDL_FOUND)
SET( LIBDL_LIBRARIES ${LIBDL_LIBRARY} )
ELSE(LIBDL_FOUND)
SET( LIBDL_LIBRARIES )
ENDIF(LIBDL_FOUND)

MARK_AS_ADVANCED( LIBDL_LIBRARY LIBDL_INCLUDE_DIR )