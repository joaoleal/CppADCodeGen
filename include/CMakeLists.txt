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
#
# ----------------------------------------------------------------------------
# files to be installed
# ----------------------------------------------------------------------------
INSTALL( DIRECTORY "cppadcg/"
         DESTINATION ${install_include_location}
         FILES_MATCHING PATTERN "*.hpp")
