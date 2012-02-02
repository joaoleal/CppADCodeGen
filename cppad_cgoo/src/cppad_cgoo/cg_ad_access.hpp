#ifndef CPPAD_CG_ACCESS_INCLUDED
#define	CPPAD_CG_ACCESS_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    template<class Base>
    class ADBaseAccess {
    protected:
        // static methods that allow access to the private data of AD<Base>

        inline static const Base& getADBaseValue(const AD<Base>& ad) {
            return ad.value_;
        }
    };

    template<class Base>
    class ADAccess : public ADBaseAccess<CG<Base> > {
    private:

        //Friends
        template<class T>
        friend AD<CG<T> > CondExpOp(enum CompareOp cop,
                                    const AD<CG<T> > &left,
                                    const AD<CG<T> > &right,
                                    const AD<CG<T> > &trueCase,
                                    const AD<CG<T> > &falseCase);

        // static methods that allow access to the private data of AD<Base>

        inline static const CG<Base>& getValue(const AD<CG<Base> >& ad) {
            return ADBaseAccess<CG<Base> >::getADBaseValue(ad);
        }
    };

}

#endif

