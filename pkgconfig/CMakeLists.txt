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
# Build the cppad pkg-config information files.
#
CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/cppadcg.pc.in             
	${CMAKE_CURRENT_BINARY_DIR}/cppadcg.pc 
)

CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/cppadcg-uninstalled.pc.in 
	${CMAKE_CURRENT_BINARY_DIR}/cppadcg-uninstalled.pc
)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/cppadcg.pc
	DESTINATION ${install_library_pkg_location} )
