#ifndef CPPAD_TESTEXCEPTION_HPP
#define CPPAD_TESTEXCEPTION_HPP
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

#include <string>
#include <stdexcept>

namespace CppAD {

    class TestException : public std::exception {
    protected:
        std::string _message;

    public:
        inline TestException(const std::string& message) throw () :
            _message(message) {
        }

        inline const char* what() const throw () {
            return _message.c_str();
        }

        inline virtual ~TestException() throw () {
        }

    private:
        TestException() throw (); // not implemented
    };

}

#endif

