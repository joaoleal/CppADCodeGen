#ifndef CPPAD_CG_CPPADCG_ASSERT_INCLUDED
#define	CPPAD_CG_CPPADCG_ASSERT_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#define CPPADCG_ASSERT_KNOWN(exp, msg)          \
{	if( ! ( exp ) )                         \
	CppAD::ErrorHandler::Call(              \
		true       ,                    \
		__LINE__   ,                    \
 		__FILE__   ,                    \
		#exp       ,                    \
		msg        );                   \
}

#endif
