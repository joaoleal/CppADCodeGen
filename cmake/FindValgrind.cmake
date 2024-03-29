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
# - Try to find valgrind
# Once done this will define
#  Valgrind_FOUND       - System has valgrind
#  Valgrind_EXECUTABLE  - The executable location

FIND_PROGRAM(Valgrind_EXECUTABLE NAME valgrind
             HINTS $ENV{Valgrind_HOME}
                   "/usr/bin/"
                   "/usr/local/bin")

# handle the QUIETLY and REQUIRED arguments and set Valgrind_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Valgrind DEFAULT_MSG Valgrind_EXECUTABLE)

MARK_AS_ADVANCED(Valgrind_EXECUTABLE)