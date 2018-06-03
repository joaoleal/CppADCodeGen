# --------------------------------------------------------------------------
#  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
#    Copyright (C) 2016 Ciengis
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
ADD_SUBDIRECTORY(dynamiclib)

IF(PDFLATEX_COMPILER)
    ADD_SUBDIRECTORY(lang/latex)
ENDIF()

IF(GRAPHVIZ_DOT_PATH)
    ADD_SUBDIRECTORY(lang/dot)
ENDIF()

ADD_SUBDIRECTORY(lang/mathml)

IF ( ENABLE_THREAD_POOL_TESTS )
    ADD_SUBDIRECTORY(threadpool)
ENDIF()

IF(LLVM_FOUND AND "${LLVM_VERSION_MAJOR}.${LLVM_VERSION_MINOR}" MATCHES "^(3.2|3.3|3.4|3.6|3.8|4.0|5.0|6.0)$")
    IF((("${LLVM_VERSION_MAJOR}.${LLVM_VERSION_MINOR}" MATCHES "^(3.2|3.6|3.8|4.0|5.0|6.0)$") AND CLANG_FOUND)
            OR
    ("${LLVM_VERSION_MAJOR}.${LLVM_VERSION_MINOR}" MATCHES "^(3.3|3.4)$")
            )
        ADD_SUBDIRECTORY(llvm)
    ENDIF()
ENDIF()
