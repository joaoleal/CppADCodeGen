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

