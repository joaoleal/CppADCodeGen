# ----------------------------------------------------------------------------
# CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
#   Copyright (C) 2012 Ciengis
#
# CppADCodeGen is distributed under multiple licenses:
#
#  - Common Public License Version 1.0 (CPL1), and
#  - GNU General Public License Version 2 (GPL2).
#
# CPL1 terms and conditions can be found in the file epl-v10.txt, while
# terms and conditions for the GPL2 can be found in the file gpl2.txt.
# ----------------------------------------------------------------------------
#
# - Try to find CppAD
# Once done this will define
#  CPPAD_FOUND - System has CppAD
#  CPPAD_INCLUDE_DIRS - The CppAD include directories
#  CPPAD_LIBRARY_DIRS - The library directories needed to use CppAD
#  CPPAD_LIBRARIES    - The libraries needed to use CppAD

IF (CPPAD_INCLUDES AND CPPAD_LIBRARIES)
  SET(CPPAD_FIND_QUIETLY TRUE)
ENDIF (CPPAD_INCLUDES AND CPPAD_LIBRARIES)


FIND_PACKAGE(PkgConfig)
IF( PKG_CONFIG_FOUND )
  pkg_check_modules(CPPAD QUIET cppad)
ENDIF( PKG_CONFIG_FOUND )


IF( CPPAD_FOUND )
  IF(NOT CPPAD_FIND_QUIETLY)
    MESSAGE(STATUS "package cppad found")
  ENDIF()
ELSE( CPPAD_FOUND )
  FIND_PATH(CPPAD_INCLUDE_DIRS NAMES cppad/cppad.hpp
            HINTS  $ENV{CPPAD_HOME}
                   "/usr/include" )
           
  FIND_LIBRARY(CPPAD_IPOPT_LIBRARY 
                cppad_ipopt
                HINTS "$ENV{CPPAD_HOME}/lib"
                      "/usr/lib" )

  SET(CPPAD_INCLUDE_DIRS ${CPPAD_INCLUDE_DIR})
  SET(CPPAD_LIBRARIES ${CPPAD_IPOPT_LIBRARY})

  INCLUDE(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set LIBIPOPT_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(CPPAD  DEFAULT_MSG
                                    CPPAD_INCLUDE_DIRS)

  MARK_AS_ADVANCED(CPPAD_INCLUDE_DIRS CPPAD_LIBRARIES)
  
  IF( CPPAD_FOUND AND NOT CPPAD_FIND_QUIETLY )
    MESSAGE(STATUS "package CppAD found")
  ENDIF()
ENDIF( CPPAD_FOUND )

