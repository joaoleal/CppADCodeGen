#ifndef CPPAD_CGEXCEPTION_INCLUDED
#define	CPPAD_CGEXCEPTION_INCLUDED
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

    /**
     * Source code generation exception
     * 
     * \author Joao Leal
     */
    class CGException : public std::exception {
    protected:
        std::string _message;

    public:

        CGException(const std::string& message) throw () {
            _message = message;
        }

        const char* what() const throw () {
            return _message.c_str();
        }

        virtual ~CGException() throw () {
        }

    protected:
        CGException() throw ();
    };

}

#endif

