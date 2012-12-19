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

ADD_EXECUTABLE(source_generation EXCLUDE_FROM_ALL
	source_generation.cpp)
SET(EXAMPLES_COMMAND "COMMAND source_generation")
ADD_CUSTOM_TARGET(example_source_generation
                  COMMAND source_generation)

IF( UNIX )
	ADD_EXECUTABLE(dynamic_linux EXCLUDE_FROM_ALL
		dynamic_linux.cpp)
	SET(EXAMPLES_COMMAND "${EXAMPLES_COMMAND} COMMAND dynamic_linux")

	TARGET_LINK_LIBRARIES(dynamic_linux ${DL_LIBRARIES})
	ADD_CUSTOM_TARGET(example_dynamic_linux
		 COMMAND dynamic_linux)
ENDIF( UNIX )


ADD_CUSTOM_TARGET(examples)
                  
ADD_DEPENDENCIES(examples example_source_generation)
IF( UNIX )
	ADD_DEPENDENCIES(examples example_dynamic_linux)
ENDIF()
