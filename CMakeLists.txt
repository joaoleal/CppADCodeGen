# -----------------------------------------------------------------------------
# CppADCodeGen: C++ Algorithmic Differentiation with source code generation
#
# Copyright (C) 2012 Ciengis
#
# CppAD is distributed under multiple licenses. This distribution is under
# the terms of the 
#                     Eclipse Public License Version 1.0.
#
# A copy of this license is included in the COPYING file of this distribution.
# Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
# -----------------------------------------------------------------------------
# Some constants
# ----------------------------------------------------------------------------
# Set the minimum required version of cmake for this project.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(cppadcg CXX)

# cppadcg_version is used by set_version.sh to get the version number.
SET(cppadcg_version      "20121215" )
SET(cppadcg_url          "https://github.com/joaoleal/CppADCodeGen" )
SET(cppadcg_description  "Differentiation of C++ Algorithms with Source Code Generation" )

# some additional cmake modules
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

FIND_PACKAGE(CppAD REQUIRED)

FIND_PACKAGE(Eigen3 REQUIRED)

IF( UNIX )
	FIND_PACKAGE(DL REQUIRED)
ENDIF()

INCLUDE_DIRECTORIES( ${CPPAD_INCLUDE_DIR} "${CMAKE_SOURCE_DIR}/include" ${EIGEN3_INCLUDE_DIR} ${DL_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
#         General installation information
# ----------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX "/usr" CACHE PATH "Default install path prefix")
set(install_export_location "${CMAKE_INSTALL_PREFIX}/share/cppadcg")
set(install_include_location "${CMAKE_INSTALL_PREFIX}/include/cppadcg")
set(install_library_pkg_location "${install_library_location}/pkgconfig")


# ----------------------------------------------------------------------------
# Define some optional compiler flags
# ----------------------------------------------------------------------------
IF(CMAKE_COMPILER_IS_GNUCXX)
  
  SET(warning_flags "-Wall") # enable as many warnings as possible
  #SET(warning_flags "${warning_flags} -Wextra") 
  #SET(warning_flags "${warning_flags} -Wshadow") # definition shadows another
  SET(warning_flags "${warning_flags} -Wno-system-headers") # do not issue warnings for system headers
  #SET(warning_flags "${warning_flags} -Woverloaded-virtual") # overloaded virtual function has a different signature

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warning_flags}")  
 
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -ggdb -D_FORTIFY_SOURCE=2")

  # NDEBUG required in order to make CppAD significantly faster
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")
ENDIF()

# ----------------------------------------------------------------------------
# Currently building tests as normal executables
# ----------------------------------------------------------------------------
ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(example EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(test EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(pkgconfig)

