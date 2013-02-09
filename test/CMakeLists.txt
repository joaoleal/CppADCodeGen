# --------------------------------------------------------------------------
#  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
#    Copyright (C) 2012 Ciengis
#
#  CppADCodeGen is distributed under multiple licenses:
#
#   - Common Public License Version 1.0 (CPL1), and
#   - GNU General Public License Version 2 (GPL2).
#
# CPL1 terms and conditions can be found in the file "epl-v10.txt", while
# terms and conditions for the GPL2 can be found in the file "gpl2.txt".
# ----------------------------------------------------------------------------
#
# Author: Joao Leal
#
# ----------------------------------------------------------------------------
SET(CMAKE_BUILD_TYPE DEBUG)

################################################################################
# download and compile GTest
################################################################################
INCLUDE(ExternalProject)

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/test/googletest)

# Add gtest
ExternalProject_Add(
    googletest
    SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk/
    TIMEOUT 10
    # Force separate output paths for debug and release builds to allow easy
    # identification of correct lib in subsequent TARGET_LINK_LIBRARIES commands
    CMAKE_ARGS -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
               -Dgtest_force_shared_crt=ON
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
SET(gtest_source_dir ${source_dir})

# Specify MainTest's link libraries
ExternalProject_Get_Property(googletest binary_dir)
SET(GTEST_LIBS_DIR ${binary_dir})

INCLUDE_DIRECTORIES("${gtest_source_dir}/include")

################################################################################
# add gtest dependency
################################################################################
FUNCTION(add_gest_link_libraries test_target_name use_main)
  IF(CMAKE_COMPILER_IS_GNUCXX)
    TARGET_LINK_LIBRARIES(${test_target_name} pthread)
    TARGET_LINK_LIBRARIES(${test_target_name} ${GTEST_LIBS_DIR}/libgtest.a)
    IF(${use_main})
      TARGET_LINK_LIBRARIES(${test_target_name} ${GTEST_LIBS_DIR}/libgtest_main.a)
    ENDIF()
  ELSE()
    TARGET_LINK_LIBRARIES(${test_target_name}
                          debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
                          optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES})
    IF(${use_main})
      TARGET_LINK_LIBRARIES(${test_target_name}
                            debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
                            optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES})
    ENDIF()
  ENDIF()
  
  # Create dependency of test on googletest
  ADD_DEPENDENCIES(${test_target_name} googletest)
ENDFUNCTION()

################################################################################
# add test function
################################################################################
FUNCTION(add_cppadcg_test source_file)

  GET_FILENAME_COMPONENT(test_target_name "${source_file}" NAME_WE)

  ADD_EXECUTABLE("${test_target_name}"
                 # sources:
                 "${source_file}"
                 ${ARGN})

  # Specify MainTest's link libraries
  add_gest_link_libraries(${test_target_name} ON)

  IF( UNIX )
    TARGET_LINK_LIBRARIES(${test_target_name} ${DL_LIBRARIES})
  ENDIF()

  ADD_CUSTOM_COMMAND(TARGET ${test_target_name}
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} ARGS -E remove_directory "tmp" # make sure tmp is empty
                     COMMAND ${CMAKE_COMMAND} ARGS -E make_directory "tmp"
                     COMMENT "Creating tmp folder")

  ADD_TEST(NAME ${test_target_name}
           COMMAND ${test_target_name})
  IF(VALGRIND_FOUND)
    ADD_TEST(NAME memcheck_${test_target_name}
             COMMAND ${GTEST_MEMORYCHECK_COMMAND} 
                     --read-var-info=yes
                     --leak-check=full
                     --track-origins=yes
                     --leak-resolution=high
                     --quiet
                     --error-exitcode=1
                     --log-file=${CMAKE_CURRENT_BINARY_DIR}/valgrind.log.${test_target_name}
                     ./${test_target_name}
             COMMENT "Running valgrind on test \"${test_target_name}\"")
  ENDIF()
ENDFUNCTION()
################################################################################
# use valgrind if possible
################################################################################
FIND_PACKAGE(Valgrind)

IF(VALGRIND_FOUND)
  SET(GTEST_MEMORYCHECK_COMMAND "${VALGRIND_EXECUTABLE}")
  MESSAGE(STATUS "Valgrind found")
ELSE()
  MESSAGE(WARNING "Valgrind not found")
ENDIF()

################################################################################
# add tests
################################################################################
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

add_cppadcg_test(abs.cpp)
add_cppadcg_test(acos.cpp)
add_cppadcg_test(add.cpp)
add_cppadcg_test(asin.cpp)
add_cppadcg_test(assign.cpp)
add_cppadcg_test(atan.cpp)
add_cppadcg_test(atan_2.cpp)
add_cppadcg_test(compare_change.cpp)
add_cppadcg_test(cond_exp.cpp)
add_cppadcg_test(cos.cpp)
add_cppadcg_test(cosh.cpp)
add_cppadcg_test(div.cpp)
add_cppadcg_test(div_mul.cpp)
add_cppadcg_test(exp.cpp)
add_cppadcg_test(inputstream.cpp)
add_cppadcg_test(log.cpp)
add_cppadcg_test(log_10.cpp)
add_cppadcg_test(mul.cpp)
add_cppadcg_test(parameter.cpp)
add_cppadcg_test(pow.cpp)
add_cppadcg_test(sin.cpp)
add_cppadcg_test(sub.cpp)
add_cppadcg_test(tan.cpp)
add_cppadcg_test(unary.cpp)
IF( UNIX )
    add_cppadcg_test(dynamic.cpp)
ENDIF()

ADD_SUBDIRECTORY(solve)
ADD_SUBDIRECTORY(dae_index_reduction)

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/valgrind.sh"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")