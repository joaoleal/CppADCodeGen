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

ADD_EXECUTABLE("speed_plugflow"
               # sources:
               "job_speed_listener.cpp"
               "speed_plugflow.cpp")

IF( UNIX )
  TARGET_LINK_LIBRARIES(speed_plugflow ${DL_LIBRARIES})
ENDIF()

ADD_EXECUTABLE("speed_collocation"
               # sources:
               "job_speed_listener.cpp"
               "speed_collocation.cpp")

IF( UNIX )
  TARGET_LINK_LIBRARIES(speed_collocation ${DL_LIBRARIES})
ENDIF()