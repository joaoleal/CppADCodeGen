# ----------------------------------------------------------------------------
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
# - Try to find Clang
# Once done this will define
# Clang_FOUND        - True if Clang found.
# Clang_INCLUDE_DIRS - where to find Clang include files
# Clang_LIBS         - list of clang libs

UNSET(Clang_INCLUDE_DIRS CACHE)
UNSET(Clang_LIBS CACHE)

IF(NOT LLVM_INCLUDE_DIRS OR NOT LLVM_LIBRARY_DIRS)
  SET(Clang_FOUND FALSE)
  MESSAGE(WARNING "No LLVM and Clang support requires LLVM")
  RETURN()
ENDIF()

MACRO(findClangStaticLib _libname_)
  UNSET(Clang_${_libname_}_LIB CACHE)
  IF(LLVM_LIBRARY_DIRS)
      FIND_LIBRARY(Clang_${_libname_}_LIB ${_libname_} PATH ${LLVM_LIBRARY_DIRS})
  ELSE()
      FIND_LIBRARY(Clang_${_libname_}_LIB ${_libname_} ${LLVM_LIBRARY_DIRS} ${Clang_LIBRARY_DIRS})
  ENDIF()


  MARK_AS_ADVANCED(Clang_${_libname_}_LIB)
  IF(Clang_${_libname_}_LIB)
     SET(Clang_LIBS ${Clang_LIBS} ${Clang_${_libname_}_LIB})
  ENDIF()
ENDMACRO()

findClangStaticLib(clang NAMES clang libclang) # LibClang: high-level C interface

# Clang shared library provides just the limited C interface, so it
# can not be used.  We look for the static libraries.
SET(Clang_LIBNAMES clangCodeGen clangFrontendTool clangFrontend clangDriver clangSerialization clangTooling
                   clangParse clangSema clangChecker clangRewrite clangRewriteFrontend
                   clangStaticAnalyzerFrontend clangStaticAnalyzerCheckers
                   clangStaticAnalyzerCore clangAnalysis clangARCMigrate clangEdit clangAST clangASTMatchers clangLex
                   clangBasic clangRewriteCore)

foreach(LIBNAME ${Clang_LIBNAMES})
    findClangStaticLib(${LIBNAME})
endforeach()

UNSET(Clang_INCLUDE_DIRS CACHE)
FIND_PATH(Clang_INCLUDE_DIRS clang/Basic/Version.h HINTS ${LLVM_INCLUDE_DIRS})

INCLUDE(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set Clang_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Clang  DEFAULT_MSG
                                  Clang_INCLUDE_DIRS
                                  Clang_LIBS)

MARK_AS_ADVANCED(Clang_INCLUDE_DIRS Clang_LIBS)
