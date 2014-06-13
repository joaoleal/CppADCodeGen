# --------------------------------------------------------------------------
#  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
#    Copyright (C) 2012 Ciengis
#
#  CppADCodeGen is distributed under multiple licenses:
#
#   - Eclipse Public License Version 1.0 (EPL1), and
#   - GNU General Public License Version 3 (GPL3).
#
# EPL1 terms and conditions can be found in the file "epl-v10.txt", while
# terms and conditions for the GPL3 can be found in the file "gpl3.txt".
# ----------------------------------------------------------------------------
#
# Author: Joao Leal
#
# ----------------------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(cppadcg CXX C)

SET(cppadcg_version      "2.0.0" )
SET(cppadcg_url          "https://github.com/joaoleal/CppADCodeGen" )
SET(cppadcg_description  "A C++ Algorithmic Differentiation Package with Source Code Generation" )

# some additional cmake modules
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

FIND_PACKAGE(Threads REQUIRED)

FIND_PACKAGE(CppAD REQUIRED)

FIND_PACKAGE(ADOLC)

FIND_PACKAGE(Eigen3)

IF(UNIX)
    FIND_PACKAGE(DL REQUIRED)
ENDIF()

FIND_PACKAGE(LLVM)
FIND_PACKAGE(Clang)

IF(LLVM_FOUND)
    ADD_DEFINITIONS(-DLLVM_VERSION_MAJOR=${LLVM_VERSION_MAJOR}
                    -DLLVM_VERSION_MINOR=${LLVM_VERSION_MINOR})
ENDIF()

###
# cppcheck tests
###
INCLUDE(cmake/cppcheck.cmake)

################################################################################
#                              Includes
################################################################################
INCLUDE_DIRECTORIES(${CPPAD_INCLUDE_DIR} "${CMAKE_SOURCE_DIR}/include" ${DL_INCLUDE_DIRS})

IF(EIGEN3_FOUND)
    INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
ENDIF()

################################################################################
#         General installation information
################################################################################
IF(UNIX)
    SET(CMAKE_INSTALL_PREFIX "/usr" CACHE PATH "Default install path prefix")
ENDIF()

SET(install_cppadcg_include_location "${CMAKE_INSTALL_PREFIX}/include/cppad")
SET(install_cppad_include_location "${CMAKE_INSTALL_PREFIX}/include/cppad")
SET(install_library_pkg_location "${CMAKE_INSTALL_PREFIX}/share/pkgconfig")


################################################################################
# Define some optional compiler flags
################################################################################
IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(warning_flags "-Wall") # enable as many warnings as possible
    #SET(warning_flags "${warning_flags} -Wextra") 
    #SET(warning_flags "${warning_flags} -Wshadow") # definition shadows another
    SET(warning_flags "${warning_flags} -Wno-system-headers") # do not issue warnings for system headers
    #SET(warning_flags "${warning_flags} -Woverloaded-virtual") # overloaded virtual function has a different signature

    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${warning_flags}")  

    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -ggdb -D_FORTIFY_SOURCE=2")

    # NDEBUG required in order to make CppAD significantly faster
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")

    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -ggdb -DNDEBUG -D_FORTIFY_SOURCE=2")
ENDIF()

################################################################################
# Currently building tests as normal executables
################################################################################
ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(example EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(speed EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(pkgconfig)

# enable testing must be called here in order for the target test to be created
ENABLE_TESTING()
OPTION(ENABLE_TEST_CPPCHECKS "Use CppCheck to verify source code in tests (very lengthy task!)" OFF)
ADD_SUBDIRECTORY(test EXCLUDE_FROM_ALL)

################################################################################
# Generate the API documentation with Doxygen
################################################################################
FIND_PACKAGE(Doxygen)

IF(DOXYGEN_FOUND)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doxyfile @ONLY)
    ADD_CUSTOM_TARGET(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating API documentation with Doxygen" VERBATIM)
ENDIF()
