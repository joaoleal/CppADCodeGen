/* 
 * File:   main.cpp
 * Author: joao
 *
 * Created on 21 de Janeiro de 2012, 16:57
 */

#include <cstdlib>
#include <vector>

#include <cg.hpp>

using namespace std;
using namespace CppAD;

/*
 * 
 */
int main(int argc, char** argv) {

    CodeHandler<double> handler(std::cout);

    std::vector<CG<double> > v(6);
    handler.makeVariables(v);

    std::vector<AD<CG<double> > > V(v.size());
    for (size_t i = 0; i < V.size(); i++) {
        V[i] = v[i];
    }

    CppAD::Independent(V);

    // dependent variable values
    std::vector<AD<CG<double> > > z(3);

    z[0] = 0.0 + v[1] / 1. + v[2] + 2 - (v[1] + v[2]) + 5.0 / (v[3] + v[4]) * v[5] * v[4] * 1.0;
    z[1] = v[1];
    z[0] = sin(v[0]);

    //    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CG<double> > f(V, z);


    cout << V[0];
    
    return EXIT_SUCCESS;
}

