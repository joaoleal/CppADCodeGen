#ifndef CPPADCGOO_TEST_ADD_INCLUDED
#define	CPPADCGOO_TEST_ADD_INCLUDED

#include <assert.h>

template<class T>
CppAD::ADFun<T>* AddFunc(const std::vector<CppAD::AD<T> >& u) {
    using namespace CppAD;
    using namespace std;

    assert(u.size() == 2);

    size_t s = 0;
    size_t t = 1;

    // dependent variable vector and indices
    std::vector< AD<T> > Z(3);
    size_t x = 0;
    size_t y = 1;
    size_t z = 2;

    // dependent variable values
    Z[x] = u[s] + u[t]; // AD<double> + AD<double>
    Z[y] = Z[x] + 1.; // AD<double> + double
    Z[z] = 1. + Z[y]; // double + AD<double> 

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (u, Z);
}

#endif

