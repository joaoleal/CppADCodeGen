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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${LLVM_INCLUDE_DIRS} ${DL_INCLUDE_DIRS})

create_temp_folder_target()

LINK_DIRECTORIES(${LLVM_LIBRARY_DIRS})
ADD_DEFINITIONS(${LLVM_CFLAGS})

add_cppadcg_test(llvm.cpp)

IF(${LLVM_VERSION_MAJOR} STREQUAL 3 AND ${LLVM_VERSION_MINOR} STREQUAL 2)
  TARGET_LINK_LIBRARIES(llvm
                        ${CLANG_LIBS}
                        ${LLVM_MODULE_LIBS}
                        ${LLVM_LDFLAGS})
ELSE()
  TARGET_LINK_LIBRARIES(llvm
                        ${LLVM_MODULE_LIBS}
                        ${LLVM_LDFLAGS})
ENDIF()
