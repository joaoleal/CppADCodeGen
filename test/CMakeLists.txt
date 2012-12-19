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
#
# Specifies build type for this directory. Possible values are 
# empty, Debug, Release, RelWithDebInfo and MinSizeRel
SET(CMAKE_BUILD_TYPE DEBUG)

IF( UNIX )
	SET(UNIX_TESTS dynamic.cpp)
ENDIF()

ADD_EXECUTABLE(testCGOO EXCLUDE_FROM_ALL
	#sources:
	gcc_load_dynamic.hpp
	TestException.hpp
	test_code_gen.cpp
	abs.cpp
	acos.cpp
	add.cpp
	asin.cpp
	assign.cpp
	atan.cpp
	atan_2.cpp
	compare_change.cpp
	cond_exp.cpp
	cos.cpp
	cosh.cpp
	div.cpp
	exp.cpp
	log.cpp
	log_10.cpp
	mul.cpp
	parameter.cpp
	pow.cpp
	sin.cpp
	sub.cpp
	tan.cpp
	unary.cpp
	gcc_load_dynamic.cpp
	${UNIX_TESTS}
	solve/solve_add.cpp
	solve/solve_cosh.cpp
	solve/solve_div.cpp
	solve/solve_exp.cpp
	solve/solve_log.cpp
	solve/solve_log_10.cpp
	solve/solve_mul.cpp
	solve/solve_pow.cpp
	solve/solve_sinh.cpp
	solve/solve_sqrt.cpp
	solve/solve_sub.cpp
	solve/solve_tanh.cpp
	solve/test_solve.cpp
	solve/solve_unary.cpp
	dae_index_reduction/dummy_derivative.cpp
	dae_index_reduction/pantelides.cpp
	TestException.cpp
)

IF( UNIX )
	TARGET_LINK_LIBRARIES(testCGOO ${DL_LIBRARIES})
ENDIF()

ADD_CUSTOM_TARGET(test
                  COMMAND ${CMAKE_COMMAND} -E remove_directory "tmp" # make sure tmp is empty
                  COMMAND ${CMAKE_COMMAND} -E make_directory "tmp"
                  COMMAND testCGOO
                  COMMAND ${CMAKE_COMMAND} -E remove_directory "tmp")
                  
ADD_DEPENDENCIES(test testCGOO)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_TEST(NAME testCGOO
         COMMAND testCGOO)