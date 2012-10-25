/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include "TestException.hpp"

using namespace CppAD;

TestException::TestException() throw () {
}

TestException::TestException(const std::string& message) throw () {
    _message = message;
}

const char* TestException::what() const throw () {
    return _message.c_str();
}

TestException::~TestException() throw () {
}

