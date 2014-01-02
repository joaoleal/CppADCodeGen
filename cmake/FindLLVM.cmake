# ----------------------------------------------------------------------------
#  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
#    Copyright (C) 2013 Ciengis
#
#  CppADCodeGen is distributed under multiple licenses:
#
#   - Eclipse Public License Version 1.0 (EPL1), and
#   - GNU General Public License Version 3 (GPL3).
#
#  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
#  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
# ----------------------------------------------------------------------------
#
# - Try to find LLVM
#
# It defines the following variables
#  LLVM_FOUND        - True if llvm found.
#  LLVM_INCLUDE_DIRS - where to find llvm include files
#  LLVM_LIBRARY_DIRS - where to find llvm libs
#  LLVM_CFLAGS       - llvm compiler flags
#  LLVM_LDFLAGS      - llvm linker flags
#  LLVM_MODULE_LIBS  - list of llvm libs for working with modules.

FIND_PROGRAM(LLVM_CONFIG llvm-config)

IF(LLVM_CONFIG)
  MESSAGE(STATUS "llvm-config found at: ${LLVM_CONFIG}")
ELSE()
  MESSAGE(FATAL_ERROR "Could NOT find llvm-config")
ENDIF()

EXECUTE_PROCESS(COMMAND ${LLVM_CONFIG} --version
                OUTPUT_VARIABLE LLVM_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)

STRING(REGEX REPLACE "^([0-9]+)\\.([0-9]+).*" "\\1"
       LLVM_VERSION_MAJOR
       "${LLVM_VERSION}")

STRING(REGEX REPLACE "^([0-9]+)\\.([0-9]+).*" "\\2"
       LLVM_VERSION_MINOR
       "${LLVM_VERSION}")

EXECUTE_PROCESS(COMMAND ${LLVM_CONFIG} --includedir
                OUTPUT_VARIABLE LLVM_INCLUDE_DIRS
                OUTPUT_STRIP_TRAILING_WHITESPACE)

EXECUTE_PROCESS(COMMAND ${LLVM_CONFIG} --libdir
                OUTPUT_VARIABLE LLVM_LIBRARY_DIRS
                OUTPUT_STRIP_TRAILING_WHITESPACE)

EXECUTE_PROCESS(COMMAND ${LLVM_CONFIG} --cppflags
                OUTPUT_VARIABLE LLVM_CFLAGS
                OUTPUT_STRIP_TRAILING_WHITESPACE)

IF(LLVM_CFLAGS MATCHES "\\-DNDEBUG")
    SET(LLVM_WITH_NDEBUG TRUE)
ELSE()
    SET(LLVM_WITH_NDEBUG FALSE)
ENDIF()

FIND_LIBRARY(LLVM_MODULE_LIBS LLVM-${LLVM_VERSION_MAJOR}.${LLVM_VERSION_MINOR} ${LLVM_LIBRARY_DIRS})
IF(NOT LLVM_MODULE_LIBS)
  EXECUTE_PROCESS(COMMAND ${LLVM_CONFIG} --libs
                  OUTPUT_VARIABLE LLVM_MODULE_LIBS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
ENDIF()

EXECUTE_PROCESS(COMMAND ${LLVM_CONFIG} --ldflags
                OUTPUT_VARIABLE LLVM_LDFLAGS
                OUTPUT_STRIP_TRAILING_WHITESPACE)

INCLUDE(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LLVM_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LLVM  DEFAULT_MSG
                                  LLVM_INCLUDE_DIRS
                                  LLVM_LIBRARY_DIRS
                                  LLVM_CFLAGS
                                  LLVM_LDFLAGS
                                  LLVM_MODULE_LIBS)

MARK_AS_ADVANCED(LLVM_INCLUDE_DIRS
                 LLVM_LIBRARY_DIRS
                 LLVM_CFLAGS
                 LLVM_LDFLAGS
                 LLVM_MODULE_LIBS)