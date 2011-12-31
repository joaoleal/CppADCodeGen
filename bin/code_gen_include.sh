#! /bin/bash -e
# $Id: search.sh 2082 2011-08-31 17:50:58Z bradbell $
# -----------------------------------------------------------------------------
# CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell
#
# CppAD is distributed under multiple licenses. This distribution is under
# the terms of the
#                     Common Public License Version 1.0.
#
# A copy of this license is included in the COPYING file of this distribution.
# Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
# -----------------------------------------------------------------------------
if [ ! -e "bin/code_gen_include.sh" ]
then
	echo "bin/code_gen_include.sh: must be executed from its parent directory"
	exit 1
fi
#
for file in cppad/local/code_gen/*.hpp
do
	sed < $file > junk.$$ \
		-e 's|^\(# *include *\)"\([^"]*\)"|\1 <cppad/local/code_gen/\2>|' 
	echo $file
	if  ! diff $file junk.$$
	then
		mv junk.$$ $file
	fi
done
rm junk.$$
