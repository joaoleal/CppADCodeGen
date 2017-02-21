/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */
#include <vector>
#include <valarray>

#include <gtest/gtest.h>

#include <cppad/cg.hpp>
#include <cppad/utility/memory_leak.hpp>
#include <cppad/utility/near_equal.hpp>

#ifndef CPPAD_CG_CPPADCGTEST_HPP
#define	CPPAD_CG_CPPADCGTEST_HPP

namespace CppAD {
namespace cg {

class CppADCGTest : public ::testing::Test {
protected:
    typedef double Base;
    typedef CppAD::cg::CG<Base> CGD;
    typedef CppAD::AD<CGD> ADCGD;
    bool verbose_;
    bool printValues_;
    bool memory_check_;
public:

    inline CppADCGTest(bool verbose = false,
                       bool printValues = false) :
        verbose_(verbose),
        printValues_(printValues),
        memory_check_(true) {
    }

    virtual void TearDown() {
        if (memory_check_) {
            ASSERT_FALSE(CppAD::memory_leak());
        }
    }

protected:

    template<class T>
    static inline ::testing::AssertionResult compareValues(const std::vector<std::vector<T> >& depCGen,
                                                           const std::vector<std::vector<T> >& dep,
                                                           T epsilonR = std::numeric_limits<T>::epsilon() * 100,
                                                           T epsilonA = std::numeric_limits<T>::epsilon() * 100) {

        assert(depCGen.size() == dep.size());

        for (size_t i = 0; i < depCGen.size(); i++) {
            ::testing::AssertionResult r = compareValues<T>(depCGen[i], dep[i], epsilonR, epsilonA);
            if (!r) {
                return ::testing::AssertionFailure() << "Comparison failed for array " << i << ":\n" << r.failure_message();
            }
        }

        return ::testing::AssertionSuccess();
    }

    template<class T>
    static inline ::testing::AssertionResult compareValues(const std::vector<T>& cgen,
                                                           const std::vector<T>& orig,
                                                           T epsilonR = std::numeric_limits<T>::epsilon() * 100,
                                                           T epsilonA = std::numeric_limits<T>::epsilon() * 100) {
        return compareValues<std::vector<T>, std::vector<T>, T>(cgen, orig, epsilonR, epsilonA);
    }

    template<class T>
    static inline ::testing::AssertionResult compareValues(const CppAD::vector<T>& cgen,
                                                           const CppAD::vector<T>& orig,
                                                           T epsilonR = std::numeric_limits<T>::epsilon() * 100,
                                                           T epsilonA = std::numeric_limits<T>::epsilon() * 100) {
        return compareValues<CppAD::vector<T>, CppAD::vector<T>, T>(cgen, orig, epsilonR, epsilonA);
    }

    template<class T>
    static inline ::testing::AssertionResult compareValues(const std::valarray<T>& cgen,
                                                           const std::valarray<T>& orig,
                                                           T epsilonR = std::numeric_limits<T>::epsilon() * 100,
                                                           T epsilonA = std::numeric_limits<T>::epsilon() * 100) {
        return compareValues<std::valarray<T>, std::valarray<T>, T>(cgen, orig, epsilonR, epsilonA);
    }

    template<class VectorT, class VectorT2, class T>
    static inline ::testing::AssertionResult compareValues(const VectorT& cgen,
                                                           const VectorT2& orig,
                                                           T epsilonR = std::numeric_limits<T>::epsilon() * 100,
                                                           T epsilonA = std::numeric_limits<T>::epsilon() * 100) {
        std::ostringstream ss;
        for (size_t i = 0; i < cgen.size(); i++) {
            ::testing::AssertionResult r = nearEqual(cgen[i], orig[i], epsilonR, epsilonA);
            if (!r) {
                ss << "Failed at element " << i << ": " << r.failure_message() << "\n";
            }
        }

        if (ss.str().empty())
            return ::testing::AssertionSuccess();
        else
            return ::testing::AssertionFailure() << ss.str();
    }

    template<class VectorBase, class T>
    static inline ::testing::AssertionResult compareValues(const VectorBase& depCGen,
                                                           const CppAD::vector<CppAD::cg::CG<T> >& dep,
                                                           T epsilonR = std::numeric_limits<T>::epsilon() * 100,
                                                           T epsilonA = std::numeric_limits<T>::epsilon() * 100) {

        std::vector<T> depd(dep.size());

        for (size_t i = 0; i < depd.size(); i++) {
            depd[i] = dep[i].getValue();
        }

        return compareValues(depCGen, &depd[0], epsilonR, epsilonA);
    }

    template<class VectorBase, class T>
    static inline ::testing::AssertionResult compareValues(const VectorBase& depCGen,
                                                           const std::vector<CppAD::cg::CG<T> >& dep,
                                                           T epsilonR = std::numeric_limits<T>::epsilon() * 100,
                                                           T epsilonA = std::numeric_limits<T>::epsilon() * 100) {

        std::vector<T> depd(dep.size());

        for (size_t i = 0; i < depd.size(); i++) {
            depd[i] = dep[i].getValue();
        }

        return compareValues(depCGen, &depd[0], epsilonR, epsilonA);
    }

    template<class VectorBool>
    static inline void compareBoolValues(const VectorBool& expected, const VectorBool& value) {
        ASSERT_EQ(expected.size(), value.size());
        for (size_t i = 0; i < expected.size(); i++) {
            ASSERT_EQ(expected[i], value[i]);
        }
    }

    template<class VectorSet>
    static inline void compareVectorSetValues(const VectorSet& expected, const VectorSet& value) {
        ASSERT_EQ(expected.size(), value.size());
        for (size_t i = 0; i < expected.size(); i++) {
            ASSERT_EQ(expected[i].size(), value[i].size());
            std::set<size_t>::const_iterator itE = expected[i].begin();
            std::set<size_t>::const_iterator itV = value[i].begin();
            for (; itE != expected[i].end(); ++itE, ++itV) {
                ASSERT_EQ(*itE, *itV);
            }
        }
    }

    template <class T>
    static inline ::testing::AssertionResult nearEqual(const T &x, const T &y,
                                                       const T &r = std::numeric_limits<T>::epsilon() * 100,
                                                       const T &a = std::numeric_limits<T>::epsilon() * 100) {

        T zero(0);
        if (r <= zero)
            return ::testing::AssertionFailure() << "Invalid relative tolerance. Must be positive!";
        if (a <= zero)
            return ::testing::AssertionFailure() << "Invalid absolute tolerance. Must be positive!";

        // check for special cases
        T infinity = T(1) / zero;

        // NaN case
        if (x != x)
            return ::testing::AssertionFailure() << "NaN value found!";
        if (y != y)
            return ::testing::AssertionFailure() << "NaN value found!";

        // infinite cases
        if (x == infinity || x == -infinity)
            return ::testing::AssertionFailure() << "Infinity value found!";
        if (y == infinity || y == -infinity)
            return ::testing::AssertionFailure() << "Infinity value found!";

        T ax = x;
        if (ax < zero)
            ax = -ax;

        T ay = y;
        if (ay < zero)
            ay = -ay;

        T ad = x - y;
        if (ad < zero)
            ad = -ad;

        if (ad > a) {
            T rel = r * (ax + ay);
            if (ad > rel) {// compare like this to avoid problems when ax and ay are zero
                return ::testing::AssertionFailure() << "Absolute error (" << ad << ") is higher than tolerance (" << a << ")"
                        "\n       and relative error (" << (ad / (ax + ay)) << ") is higher than tolerance (" << r << ")"
                        "\n       for the values " << x << " and " << y;
            }
        }

        return ::testing::AssertionSuccess();
    }
};

} // END cg namespace
} // END CppAD namespace

#endif
