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
