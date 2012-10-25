#ifndef CPPAD_TESTEXCEPTION_HPP
#define	CPPAD_TESTEXCEPTION_HPP
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <string>
#include <stdexcept>

namespace CppAD {

    class TestException : public std::exception {
    protected:
        std::string _message;

    public:
        TestException(const std::string& message) throw ();

        const char* what() const throw ();

        ~TestException() throw ();

    protected:
        TestException() throw ();
    };

}

#endif

