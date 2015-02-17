# --------------------------------------------------------------------------
#  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
#    Copyright (C) 2012 Ciengis
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
# Author: Joao Leal
#
# ----------------------------------------------------------------------------

ADD_EXECUTABLE(source_generation source_generation.cpp)
SET(EXAMPLES_COMMAND "COMMAND source_generation")
ADD_CUSTOM_TARGET(example_source_generation
                  COMMAND source_generation)

IF(UNIX)
    ADD_EXECUTABLE(dynamic_linux dynamic_linux.cpp)
    SET(EXAMPLES_COMMAND "${EXAMPLES_COMMAND} COMMAND dynamic_linux")
    TARGET_LINK_LIBRARIES(dynamic_linux ${DL_LIBRARIES})
    ADD_CUSTOM_TARGET(example_dynamic_linux
                      COMMAND dynamic_linux)

    ADD_EXECUTABLE(lagrangian lagrangian.cpp)
    SET(EXAMPLES_COMMAND "${EXAMPLES_COMMAND} COMMAND lagrangian")
    TARGET_LINK_LIBRARIES(lagrangian ${DL_LIBRARIES})
    ADD_CUSTOM_TARGET(example_lagrangian
                      COMMAND lagrangian)

    # 
    ADD_EXECUTABLE(patterns patterns.cpp)
    SET(EXAMPLES_COMMAND "${EXAMPLES_COMMAND} COMMAND patterns")
    TARGET_LINK_LIBRARIES(patterns ${DL_LIBRARIES})
    ADD_CUSTOM_TARGET(example_patterns
                      COMMAND patterns)
ENDIF(UNIX)


ADD_CUSTOM_TARGET(examples)
                  
ADD_DEPENDENCIES(examples example_source_generation)
IF(UNIX)
    ADD_DEPENDENCIES(examples example_dynamic_linux)
    ADD_DEPENDENCIES(examples example_lagrangian)
ENDIF()
