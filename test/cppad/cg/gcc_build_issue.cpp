/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Joris Vaillant
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joris Vaillant
 */
#include "CppADCGTest.hpp"

using namespace CppAD;
using namespace CppAD::cg;

// Test GCC compilation bug issue: https://github.com/joaoleal/CppADCodeGen/issues/92
// This test only have to build to pass

struct A {};
struct B {};

template <typename T>
struct Wrapper {
  virtual void foo() { 
    v.FUNC();
  }
  CppAD::cg::CGOStreamFunc<T> v;
};

int main() {
};

TEST_F(CppADCGTest, TestGCCBuildIssue) {
  Wrapper<A> a;
  Wrapper<B> b;
  a.v.FUNC();
}
