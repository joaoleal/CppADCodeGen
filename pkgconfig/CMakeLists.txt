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
# Build the cppadcg pkg-config information files.
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
