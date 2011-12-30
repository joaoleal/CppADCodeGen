#ifndef CPPAD_TESTEXCEPTION_HPP
#define	CPPAD_TESTEXCEPTION_HPP

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

