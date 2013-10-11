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

create_temp_folder_target()

add_cppadcg_test(inputstream.cpp)
add_cppadcg_test(temporary.cpp)
add_cppadcg_test(mult_sparsity_pattern.cpp)

ADD_SUBDIRECTORY(operations)
ADD_SUBDIRECTORY(models)
ADD_SUBDIRECTORY(dynamiclib)
ADD_SUBDIRECTORY(patterns)
ADD_SUBDIRECTORY(solve)
ADD_SUBDIRECTORY(dae_index_reduction)
