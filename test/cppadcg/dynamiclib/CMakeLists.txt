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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

create_temp_folder_target()

IF( UNIX )
    add_cppadcg_test(dynamic.cpp)
    add_cppadcg_test(dynamic_atomic.cpp)
    add_cppadcg_test(dynamic_atomic_2.cpp)
    add_cppadcg_test(dynamic_atomic_3.cpp)
    add_cppadcg_test(dynamic_forward_reverse.cpp)
    add_cppadcg_test(dynamic_forward_reverse_2.cpp)
ENDIF()
