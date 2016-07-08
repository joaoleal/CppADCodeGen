# --------------------------------------------------------------------------
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
# Author: Joao Leal
#
# ----------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

add_cppadcg_test(inputstream.cpp)
add_cppadcg_test(temporary.cpp)
add_cppadcg_test(mult_sparsity_pattern.cpp)

ADD_SUBDIRECTORY(extra)
ADD_SUBDIRECTORY(operations)
ADD_SUBDIRECTORY(models)
ADD_SUBDIRECTORY(dynamiclib)
ADD_SUBDIRECTORY(patterns)
ADD_SUBDIRECTORY(evaluator)
ADD_SUBDIRECTORY(solve)
ADD_SUBDIRECTORY(dae_index_reduction)

IF(PDFLATEX_COMPILER)
    ADD_SUBDIRECTORY(latex)
ENDIF()

IF(GRAPHVIZ_PATH)
    ADD_SUBDIRECTORY(dot)
ENDIF()

ADD_SUBDIRECTORY(mathml)

IF(LLVM_FOUND AND ${LLVM_VERSION_MAJOR} STREQUAL 3)
  IF(((${LLVM_VERSION_MINOR} STREQUAL 2 OR ${LLVM_VERSION_MINOR} STREQUAL 6) AND CLANG_FOUND)
     OR 
     (${LLVM_VERSION_MINOR} STREQUAL 3 OR ${LLVM_VERSION_MINOR} STREQUAL 4)
    )
    ADD_SUBDIRECTORY(llvm)
  ENDIF()
ENDIF()
