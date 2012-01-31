#ifndef CPPAD_CG_DECLARE_CG_INCLUDED
#define	CPPAD_CG_DECLARE_CG_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

// forward declarations
namespace CppAD {

    template<class Base>
    class CodeHandler;

    template<class Base>
    class CG;

    // order determining functions, see ordered.hpp
    template<class Base>
    bool GreaterThanZero(const CG<Base> &x);

    template<class Base>
    bool GreaterThanOrZero(const CG<Base> &x);

    template<class Base>
    bool LessThanZero(const CG<Base> &x);

    template<class Base>
    bool LessThanOrZero(const CG<Base> &x);

    template<class Base>
    bool abs_geq(const CG<Base>& x, const CG<Base>& y);

    // The identical property functions, see identical.hpp
    template<class Base>
    bool IdenticalZero(const CG<Base> &x) throw (CGException);

    template<class Base>
    bool IdenticalOne(const CG<Base> &x) throw (CGException);

    template<class Base>
    bool IdenticalEqualPar(const CG<Base> &x, const CG<Base> &y) throw (CGException);

    // EqualOpSeq function
    template<class Base>
    bool EqualOpSeq(const CG<Base> &u, const CG<Base> &v);

    // NearEqual function
    template<class Base>
    bool NearEqual(const CG<Base> &x, const CG<Base> &y, const Base &r, const Base &a);

    template<class Base>
    bool NearEqual(const Base &x, const CG<Base> &y, const Base &r, const Base &a);

    template<class Base>
    bool NearEqual(const CG<Base> &x, const Base &y, const Base &r, const Base &a);

    // CondExp function
    template<class Base>
    CG<Base> CondExpOp(
    enum CompareOp cop,
    const CG<Base> &left,
    const CG<Base> &right,
    const CG<Base> &trueCase,
    const CG<Base> &falseCase
    );

    //    template<class Base>
    //    bool operator<(const CG<Base> &left, const CG<Base> &right);
    //
    //    template<class Base>
    //    bool operator <=(const CG<Base> &left, const CG<Base> &right);
    //
    //    template<class Base>
    //    bool operator>(const CG<Base> &left, const CG<Base> &right);
    //
    //    template<class Base>
    //    bool operator >=(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    bool operator ==(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    bool operator !=(const CG<Base> &left, const CG<Base> &right);


    /**
     * Math functions
     */
    template<class Base>
    inline CG<Base> sign(const CG<Base> &x);

    // power function
    template<class Base>
    inline CG<Base> pow(const CG<Base> &x, const CG<Base> &y);

    template<class Base>
    inline CG<Base> abs(CG<Base> var);

    template<class Base>
    inline CG<Base> acos(CG<Base> var);

    template<class Base>
    inline CG<Base> asin(CG<Base> var);

    template<class Base>
    inline CG<Base> atan(CG<Base> var);

    template<class Base>
    inline CG<Base> cos(CG<Base> var);

    template<class Base>
    inline CG<Base> cosh(CG<Base> var);

    template<class Base>
    inline CG<Base> exp(CG<Base> var);

    template<class Base>
    inline CG<Base> log(CG<Base> var);

    template<class Base>
    inline CG<Base> sin(CG<Base> var);

    template<class Base>
    inline CG<Base> sinh(CG<Base> var);

    template<class Base>
    inline CG<Base> sqrt(CG<Base> var);

    template<class Base>
    inline CG<Base> tan(CG<Base> var);

    template<class Base>
    inline CG<Base> tanh(CG<Base> var);

}

#endif

